DROP TABLE IF EXISTS project_d2v2.d2v2_mcs_stays;
CREATE TABLE project_d2v2.d2v2_mcs_stays AS
WITH procedure_stays AS (
    SELECT DISTINCT p.stay_id, m.mcs_type
    FROM mimiciv_icu.procedureevents AS p
    INNER JOIN project_d2v2.d2v2_map_mcs AS m
      ON p.itemid = m.itemid
     AND m.linksto = 'procedureevents'
     AND m.include_main
),
chart_stays AS (
    SELECT DISTINCT c.stay_id, m.mcs_type
    FROM mimiciv_icu.chartevents AS c
    INNER JOIN project_d2v2.d2v2_map_mcs AS m
      ON c.itemid = m.itemid
     AND m.linksto = 'chartevents'
     AND m.include_main
)
SELECT
    stay_id,
    string_agg(DISTINCT mcs_type, ';' ORDER BY mcs_type) AS mcs_types,
    bool_or(source = 'procedureevents') AS has_procedure_source,
    bool_or(source = 'chartevents') AS has_chart_source
FROM (
    SELECT stay_id, mcs_type, 'procedureevents'::text AS source
    FROM procedure_stays
    UNION ALL
    SELECT stay_id, mcs_type, 'chartevents'::text AS source
    FROM chart_stays
) AS u
GROUP BY stay_id;

CREATE UNIQUE INDEX d2v2_mcs_stays_idx
    ON project_d2v2.d2v2_mcs_stays (stay_id);

-- Admission-level ICD exclusions.
DROP TABLE IF EXISTS project_d2v2.d2v2_icd_flags;
CREATE TABLE project_d2v2.d2v2_icd_flags AS
SELECT
    d.hadm_id,
    count(*) FILTER (
        WHERE m.shock_class = 'CARDIOGENIC_SHOCK'
    ) > 0 AS cardiogenic_icd,
    count(*) FILTER (
        WHERE m.shock_class = 'HEMORRHAGIC_OR_HYPOVOLEMIC_SHOCK'
    ) > 0 AS hemorrhagic_icd
FROM mimiciv_hosp.diagnoses_icd AS d
INNER JOIN project_d2v2.d2v2_map_shock_icd AS m
  ON d.icd_code = m.icd_code
 AND d.icd_version = m.icd_version
 AND m.include_main
GROUP BY d.hadm_id;

CREATE UNIQUE INDEX d2v2_icd_flags_idx
    ON project_d2v2.d2v2_icd_flags (hadm_id);

