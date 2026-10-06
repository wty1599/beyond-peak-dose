sap1_frozen_cohort AS (SELECT r.*,a.race,p.anchor_year_group,p.gender,
    (p.anchor_age+extract(year from a.admittime)::integer-p.anchor_year)::integer AS age_at_admission
    FROM requested r JOIN mimiciv_hosp.admissions a USING(hadm_id) JOIN mimiciv_hosp.patients p ON r.subject_id=p.subject_id),
sap1_nee_features AS (
SELECT c.d2v2_event_id,c.peak_nee::double precision AS episode_peak_nee,
            n.norepinephrine_equivalent_dose::double precision AS last_positive_nee,
            c.positive_hours::double precision AS positive_nee_exposure_hours,c.episode_start,c.t0 AS nee_feature_time
            FROM sap1_frozen_cohort c LEFT JOIN LATERAL (
              SELECT n.* FROM mimiciv_derived.norepinephrine_equivalent_dose n
              WHERE n.stay_id=c.stay_id AND n.norepinephrine_equivalent_dose>0 AND n.endtime>c.episode_start AND n.starttime<c.t0
              ORDER BY least(n.endtime,c.t0) DESC,greatest(n.starttime,c.episode_start) DESC,n.norepinephrine_equivalent_dose DESC LIMIT 1
            ) n ON TRUE
),
sap1_peak_nee_intervals AS (
SELECT
    c.d2v2_event_id,
    c.stay_id,
    GREATEST(n.starttime, c.episode_start) AS peak_starttime,
    LEAST(n.endtime, c.t0) AS peak_endtime
FROM sap1_frozen_cohort AS c
INNER JOIN sap1_nee_features AS f
    ON c.d2v2_event_id = f.d2v2_event_id
INNER JOIN mimiciv_derived.norepinephrine_equivalent_dose AS n
    ON c.stay_id = n.stay_id
   AND n.endtime > c.episode_start
   AND n.starttime < c.t0
WHERE n.norepinephrine_equivalent_dose > 0
  AND ABS(
      n.norepinephrine_equivalent_dose::double precision
      - f.episode_peak_nee
  ) < 1e-9
  AND LEAST(n.endtime, c.t0)
      > GREATEST(n.starttime, c.episode_start)
),
sap1_multivaso_at_peak AS (
WITH vaso_intervals AS (
    SELECT
        stay_id,
        starttime,
        endtime,
        'NOREPINEPHRINE'::text AS drug
    FROM mimiciv_derived.norepinephrine
    WHERE vaso_rate > 0
    UNION ALL
    SELECT stay_id, starttime, endtime, 'EPINEPHRINE'
    FROM mimiciv_derived.epinephrine
    WHERE vaso_rate > 0
    UNION ALL
    SELECT stay_id, starttime, endtime, 'DOPAMINE'
    FROM mimiciv_derived.dopamine
    WHERE vaso_rate > 0
    UNION ALL
    SELECT stay_id, starttime, endtime, 'PHENYLEPHRINE'
    FROM mimiciv_derived.phenylephrine
    WHERE vaso_rate > 0
    UNION ALL
    SELECT stay_id, starttime, endtime, 'VASOPRESSIN'
    FROM mimiciv_derived.vasopressin
    WHERE vaso_rate > 0
),
per_peak_interval AS (
    SELECT
        p.d2v2_event_id,
        p.peak_starttime,
        p.peak_endtime,
        COUNT(DISTINCT v.drug)::integer AS active_vasopressor_types
    FROM sap1_peak_nee_intervals AS p
    LEFT JOIN vaso_intervals AS v
        ON p.stay_id = v.stay_id
       AND v.endtime > p.peak_starttime
       AND v.starttime < p.peak_endtime
    GROUP BY
        p.d2v2_event_id,
        p.peak_starttime,
        p.peak_endtime
)
SELECT
    d2v2_event_id,
    MAX(active_vasopressor_types)::integer
        AS max_active_vasopressor_types_at_peak,
    (
        MAX(active_vasopressor_types) >= 2
    ) AS at_least_two_vasopressors_at_peak,
    MIN(peak_starttime) AS first_peak_nee_time
FROM per_peak_interval
GROUP BY d2v2_event_id
),
sap1_sofa_feature AS (
SELECT c.d2v2_event_id,s.sofa_24hours::double precision AS sofa_24h_latest,s.endtime AS sofa_window_end
            FROM sap1_frozen_cohort c LEFT JOIN LATERAL (SELECT * FROM mimiciv_derived.sofa s
            WHERE s.stay_id=c.stay_id AND s.endtime<=c.t0 ORDER BY s.endtime DESC NULLS LAST LIMIT 1) s ON TRUE
),
sap1_lactate_feature AS (
SELECT
    c.d2v2_event_id,
    lactate.valuenum::double precision AS peak_lactate_12h,
    lactate.charttime AS peak_lactate_time
FROM sap1_frozen_cohort AS c
LEFT JOIN LATERAL (
    SELECT
        l.valuenum,
        l.charttime
    FROM mimiciv_hosp.labevents AS l
    INNER JOIN project_d2v2.d2v2_map_lactate AS m
        ON l.itemid = m.itemid
       AND m.include_main
    WHERE l.hadm_id = c.hadm_id
      AND l.charttime > c.t0 - INTERVAL '12 hours'
      AND l.charttime <= c.t0
      AND l.valuenum IS NOT NULL
      AND l.valuenum >= 0
    ORDER BY l.valuenum DESC, l.charttime DESC
    LIMIT 1
) AS lactate ON TRUE
),
sap1_ventilation_feature AS (
SELECT
    c.d2v2_event_id,
    COALESCE(
        BOOL_OR(v.ventilation_status = 'InvasiveVent'),
        FALSE
    ) AS invasive_mechanical_ventilation_t0,
    CASE
        WHEN BOOL_OR(v.ventilation_status = 'InvasiveVent')
            THEN c.t0
        ELSE NULL
    END AS ventilation_feature_time
FROM sap1_frozen_cohort AS c
LEFT JOIN mimiciv_derived.ventilation AS v
    ON c.stay_id = v.stay_id
   AND v.starttime <= c.t0
   AND v.endtime > c.t0
GROUP BY c.d2v2_event_id, c.t0
),
sap1_rrt_recent AS (
SELECT
    c.d2v2_event_id,
    r.charttime,
    r.dialysis_active,
    r.dialysis_type
FROM sap1_frozen_cohort AS c
INNER JOIN mimiciv_derived.rrt AS r
    ON c.stay_id = r.stay_id
   AND r.charttime > c.t0 - INTERVAL '6 hours'
   AND r.charttime <= c.t0
WHERE r.dialysis_present = 1
),
sap1_rrt_feature AS (
WITH latest_time AS (
    SELECT
        d2v2_event_id,
        MAX(charttime) AS latest_rrt_time
    FROM sap1_rrt_recent
    GROUP BY d2v2_event_id
),
latest_state AS (
    SELECT
        r.d2v2_event_id,
        t.latest_rrt_time,
        BOOL_OR(r.dialysis_active = 1) AS dialysis_active,
        STRING_AGG(
            DISTINCT r.dialysis_type,
            ';' ORDER BY r.dialysis_type
        ) FILTER (WHERE r.dialysis_type IS NOT NULL)
            AS dialysis_types
    FROM sap1_rrt_recent AS r
    INNER JOIN latest_time AS t
        ON r.d2v2_event_id = t.d2v2_event_id
       AND r.charttime = t.latest_rrt_time
    GROUP BY r.d2v2_event_id, t.latest_rrt_time
)
SELECT
    c.d2v2_event_id,
    COALESCE(s.dialysis_active, FALSE) AS rrt_active_t0,
    s.latest_rrt_time,
    s.dialysis_types
FROM sap1_frozen_cohort AS c
LEFT JOIN latest_state AS s
    ON c.d2v2_event_id = s.d2v2_event_id
),
sap1_vital_features AS (
SELECT
    c.d2v2_event_id,
    AVG(v.mbp) FILTER (
        WHERE v.mbp IS NOT NULL
          AND v.mbp_ni IS NULL
    )::double precision AS mean_invasive_map_4h,
    COUNT(v.mbp) FILTER (
        WHERE v.mbp IS NOT NULL
          AND v.mbp_ni IS NULL
    )::integer AS invasive_map_records_4h,
    MAX(v.charttime) FILTER (
        WHERE v.mbp IS NOT NULL
          AND v.mbp_ni IS NULL
    ) AS last_invasive_map_time,
    AVG(v.mbp_ni) FILTER (
        WHERE v.mbp_ni IS NOT NULL
    )::double precision AS mean_noninvasive_map_4h,
    COUNT(v.mbp_ni) FILTER (
        WHERE v.mbp_ni IS NOT NULL
    )::integer AS noninvasive_map_records_4h,
    MAX(v.charttime) FILTER (
        WHERE v.mbp_ni IS NOT NULL
    ) AS last_noninvasive_map_time,
    AVG(v.heart_rate) FILTER (
        WHERE v.heart_rate IS NOT NULL
    )::double precision AS mean_heart_rate_4h,
    COUNT(v.heart_rate) FILTER (
        WHERE v.heart_rate IS NOT NULL
    )::integer AS heart_rate_records_4h,
    MAX(v.charttime) FILTER (
        WHERE v.heart_rate IS NOT NULL
    ) AS last_heart_rate_time
FROM sap1_frozen_cohort AS c
LEFT JOIN mimiciv_derived.vitalsign AS v
    ON c.stay_id = v.stay_id
   AND v.charttime > c.t0 - INTERVAL '4 hours'
   AND v.charttime <= c.t0
GROUP BY c.d2v2_event_id
),
sap1_urine_feature AS (
SELECT
    c.d2v2_event_id,
    SUM(u.urineoutput) FILTER (
        WHERE u.urineoutput >= 0
    )::double precision AS urine_output_24h_ml,
    COUNT(u.urineoutput) FILTER (
        WHERE u.urineoutput >= 0
    )::integer AS urine_records_24h,
    MAX(u.charttime) FILTER (
        WHERE u.urineoutput >= 0
    ) AS last_urine_time
FROM sap1_frozen_cohort AS c
LEFT JOIN mimiciv_derived.urine_output AS u
    ON c.stay_id = u.stay_id
   AND u.charttime > c.t0 - INTERVAL '24 hours'
   AND u.charttime <= c.t0
GROUP BY c.d2v2_event_id
),
sap1_feature_ledger AS (
SELECT
    c.model_row_id,
    c.d2v2_event_id,
    c.subject_id,
    c.hadm_id,
    c.stay_id,
    c.t0,
    n.episode_peak_nee,
    n.last_positive_nee,
    n.positive_nee_exposure_hours,
    n.nee_feature_time,
    m.max_active_vasopressor_types_at_peak,
    m.at_least_two_vasopressors_at_peak,
    m.first_peak_nee_time,
    s.sofa_24h_latest,
    s.sofa_window_end,
    l.peak_lactate_12h,
    l.peak_lactate_time,
    v.invasive_mechanical_ventilation_t0,
    v.ventilation_feature_time,
    r.rrt_active_t0,
    r.latest_rrt_time,
    r.dialysis_types,
    CASE
        WHEN vs.invasive_map_records_4h > 0
            THEN vs.mean_invasive_map_4h
        ELSE vs.mean_noninvasive_map_4h
    END AS mean_map_4h,
    CASE
        WHEN vs.invasive_map_records_4h > 0 THEN 'INVASIVE'
        WHEN vs.noninvasive_map_records_4h > 0 THEN 'NONINVASIVE'
        ELSE NULL
    END AS map_source,
    CASE
        WHEN vs.invasive_map_records_4h > 0 THEN 1
        WHEN vs.noninvasive_map_records_4h > 0 THEN 0
        ELSE NULL
    END::integer AS map_source_invasive,
    vs.invasive_map_records_4h,
    vs.noninvasive_map_records_4h,
    CASE
        WHEN vs.invasive_map_records_4h > 0
            THEN vs.last_invasive_map_time
        ELSE vs.last_noninvasive_map_time
    END AS map_feature_time,
    vs.mean_heart_rate_4h,
    vs.heart_rate_records_4h,
    vs.last_heart_rate_time,
    u.urine_output_24h_ml,
    u.urine_records_24h,
    u.last_urine_time
FROM sap1_frozen_cohort AS c
LEFT JOIN sap1_nee_features AS n
    ON c.d2v2_event_id = n.d2v2_event_id
LEFT JOIN sap1_multivaso_at_peak AS m
    ON c.d2v2_event_id = m.d2v2_event_id
LEFT JOIN sap1_sofa_feature AS s
    ON c.d2v2_event_id = s.d2v2_event_id
LEFT JOIN sap1_lactate_feature AS l
    ON c.d2v2_event_id = l.d2v2_event_id
LEFT JOIN sap1_ventilation_feature AS v
    ON c.d2v2_event_id = v.d2v2_event_id
LEFT JOIN sap1_rrt_feature AS r
    ON c.d2v2_event_id = r.d2v2_event_id
LEFT JOIN sap1_vital_features AS vs
    ON c.d2v2_event_id = vs.d2v2_event_id
LEFT JOIN sap1_urine_feature AS u
    ON c.d2v2_event_id = u.d2v2_event_id
),
sap1_modeling_table AS (
SELECT
    c.model_row_id,
    c.d2v2_event_id,
    c.subject_id,
    c.hadm_id,
    c.stay_id,
    c.t0,
    c.outcome_restart_72h::integer AS outcome_restart_72h,
    c.secondary_death_72h::integer AS secondary_death_72h,
    c.anchor_year_group AS validation_anchor_year_group,
    c.race AS fairness_race,
    CASE
        WHEN c.age_at_admission < 65 THEN '<65'
        ELSE '>=65'
    END AS fairness_age_group,
    c.gender AS fairness_sex,
    c.age_at_admission::double precision AS predictor_age,
    CASE
        WHEN c.gender = 'F' THEN 1
        WHEN c.gender = 'M' THEN 0
        ELSE NULL
    END::integer AS predictor_female,
    f.episode_peak_nee AS predictor_episode_peak_nee,
    f.last_positive_nee AS predictor_last_positive_nee,
    f.positive_nee_exposure_hours
        AS predictor_positive_nee_hours,
    f.at_least_two_vasopressors_at_peak::integer
        AS predictor_multivaso_peak,
    f.sofa_24h_latest AS predictor_sofa_24h,
    f.peak_lactate_12h AS predictor_peak_lactate_12h,
    f.invasive_mechanical_ventilation_t0::integer
        AS predictor_mechanical_ventilation_t0,
    f.rrt_active_t0::integer AS predictor_rrt_t0,
    f.mean_map_4h AS predictor_mean_map_4h,
    f.map_source_invasive AS predictor_map_source_invasive,
    f.mean_heart_rate_4h AS predictor_mean_heart_rate_4h,
    f.urine_output_24h_ml AS predictor_urine_output_24h_ml
FROM sap1_frozen_cohort AS c
INNER JOIN sap1_feature_ledger AS f
    ON c.model_row_id = f.model_row_id
)
SELECT m.*,l.episode_peak_nee AS ledger_episode_peak_nee,l.last_positive_nee AS ledger_last_positive_nee,l.positive_nee_exposure_hours AS ledger_positive_nee_exposure_hours,l.nee_feature_time AS ledger_nee_feature_time,l.max_active_vasopressor_types_at_peak AS ledger_max_active_vasopressor_types_at_peak,l.at_least_two_vasopressors_at_peak AS ledger_at_least_two_vasopressors_at_peak,l.first_peak_nee_time AS ledger_first_peak_nee_time,l.sofa_24h_latest AS ledger_sofa_24h_latest,l.sofa_window_end AS ledger_sofa_window_end,l.peak_lactate_12h AS ledger_peak_lactate_12h,l.peak_lactate_time AS ledger_peak_lactate_time,l.invasive_mechanical_ventilation_t0 AS ledger_invasive_mechanical_ventilation_t0,l.ventilation_feature_time AS ledger_ventilation_feature_time,l.rrt_active_t0 AS ledger_rrt_active_t0,l.latest_rrt_time AS ledger_latest_rrt_time,l.dialysis_types AS ledger_dialysis_types,l.mean_map_4h AS ledger_mean_map_4h,l.map_source AS ledger_map_source,l.map_source_invasive AS ledger_map_source_invasive,l.invasive_map_records_4h AS ledger_invasive_map_records_4h,l.noninvasive_map_records_4h AS ledger_noninvasive_map_records_4h,l.map_feature_time AS ledger_map_feature_time,l.mean_heart_rate_4h AS ledger_mean_heart_rate_4h,l.heart_rate_records_4h AS ledger_heart_rate_records_4h,l.last_heart_rate_time AS ledger_last_heart_rate_time,l.urine_output_24h_ml AS ledger_urine_output_24h_ml,l.urine_records_24h AS ledger_urine_records_24h,l.last_urine_time AS ledger_last_urine_time FROM sap1_modeling_table m JOIN sap1_feature_ledger l USING(model_row_id) ORDER BY m.model_row_id;
