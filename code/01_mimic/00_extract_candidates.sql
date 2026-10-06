\set ON_ERROR_STOP on
\pset pager off
\pset footer off
BEGIN TRANSACTION ISOLATION LEVEL REPEATABLE READ READ ONLY;
SET LOCAL statement_timeout='600s';
SET LOCAL enable_nestloop=off;
SELECT current_setting('transaction_read_only') AS read_only;
\o :candidates_file
COPY (
WITH icu AS MATERIALIZED (
 SELECT i.*,a.admittime,a.dischtime,a.deathtime,a.hospital_expire_flag,a.discharge_location,
 p.dod,p.gender,ag.age,
 m.stay_id IS NOT NULL AS mcs,
 coalesce(ic.cardiogenic_icd,false) AS cardiac_icd,
 coalesce(ic.hemorrhagic_icd,false) AS hemorrhagic_icd,
 coalesce(cu.include_cardiogenic_exclusion,false) AS cardiac_unit
 FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.admissions a USING(hadm_id)
 JOIN mimiciv_hosp.patients p ON p.subject_id=i.subject_id
 JOIN mimiciv_derived.age ag ON ag.hadm_id=i.hadm_id AND ag.age>=18
 LEFT JOIN project_d2v2.d2v2_mcs_stays m ON m.stay_id=i.stay_id
 LEFT JOIN project_d2v2.d2v2_icd_flags ic ON ic.hadm_id=i.hadm_id
 LEFT JOIN project_d2v2.d2v2_map_cardiac_unit cu ON cu.first_careunit=i.first_careunit
), pos AS MATERIALIZED (
 SELECT n.stay_id,n.starttime AS raw_starttime,n.endtime AS raw_endtime,
 greatest(n.starttime,i.intime) AS starttime,least(n.endtime,i.outtime) AS endtime,
 n.norepinephrine_equivalent_dose::double precision AS nee
 FROM mimiciv_derived.norepinephrine_equivalent_dose n JOIN icu i USING(stay_id)
 WHERE n.norepinephrine_equivalent_dose>0 AND n.endtime>i.intime AND n.starttime<i.outtime
), prior AS MATERIALIZED (
 SELECT *,max(endtime) OVER(PARTITION BY stay_id ORDER BY starttime,endtime ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS prev_end
 FROM pos WHERE endtime>starttime
), island AS MATERIALIZED (
 SELECT p.*,v.variant,
 sum(CASE WHEN prev_end IS NULL OR (v.variant='new' AND starttime>=prev_end+interval '4 hours') THEN 1 ELSE 0 END)
 OVER(PARTITION BY stay_id,v.variant ORDER BY starttime,endtime ROWS UNBOUNDED PRECEDING) AS episode
 FROM prior p CROSS JOIN (VALUES('new')) v(variant)
), episode AS MATERIALIZED (
 SELECT variant,stay_id,episode,min(starttime) AS episode_start,max(endtime) AS t0,
 max(raw_endtime) AS source_last_end,
 sum(greatest(0.0,extract(epoch FROM endtime-greatest(starttime,coalesce(prev_end,starttime)))/3600.0)) AS positive_hours,
 max(nee) AS peak_nee,count(*) AS n_nee_segments,
 md5(string_agg(raw_starttime::text||'|'||raw_endtime::text||'|'||starttime::text||'|'||endtime::text||'|'||nee::text,
 ';' ORDER BY starttime,endtime,raw_starttime,raw_endtime,nee)) AS interval_composition_signature,
 lead(min(starttime)) OVER(PARTITION BY variant,stay_id ORDER BY min(starttime)) AS next_positive
 FROM island GROUP BY variant,stay_id,episode
), base AS MATERIALIZED (
 SELECT e.*,i.subject_id,i.hadm_id,i.intime,i.outtime,i.admittime,i.dischtime,i.deathtime,
 i.hospital_expire_flag,i.discharge_location,i.dod,i.first_careunit,i.mcs,i.cardiac_icd,i.hemorrhagic_icd,i.cardiac_unit
 FROM episode e JOIN icu i USING(stay_id)
 WHERE positive_hours>=6 AND peak_nee>=0.10 AND NOT i.mcs AND NOT i.cardiac_icd AND NOT i.hemorrhagic_icd AND NOT i.cardiac_unit
), ino_src AS MATERIALIZED (
 SELECT i.stay_id,i.starttime,i.endtime FROM mimiciv_icu.inputevents i
 JOIN project_d2v2.d2v2_map_inotrope m ON i.itemid=m.itemid AND m.include_main
 WHERE EXISTS(SELECT 1 FROM base b WHERE b.stay_id=i.stay_id)
), ino AS (
 SELECT b.variant,b.stay_id,b.episode,count(i.stay_id)>0 AS inotrope
 FROM base b LEFT JOIN ino_src i ON i.stay_id=b.stay_id AND i.starttime<b.t0 AND i.endtime>b.episode_start
 GROUP BY b.variant,b.stay_id,b.episode
), prbc_src AS MATERIALIZED (
 SELECT i.stay_id,i.starttime,i.amount FROM mimiciv_icu.inputevents i
 JOIN project_d2v2.d2v2_map_prbc m ON i.itemid=m.itemid AND m.include_main
 WHERE lower(coalesce(i.amountuom,''))='ml' AND i.amount>0
 AND EXISTS(SELECT 1 FROM base b WHERE b.stay_id=i.stay_id)
), prbc AS (
 SELECT b.variant,b.stay_id,b.episode,coalesce(sum(i.amount),0) AS prbc_ml
 FROM base b LEFT JOIN prbc_src i ON i.stay_id=b.stay_id
 AND i.starttime>=b.episode_start-interval '24 hours' AND i.starttime<=least(b.episode_start+interval '24 hours',b.t0)
 GROUP BY b.variant,b.stay_id,b.episode
), lac_src AS MATERIALIZED (
 SELECT l.hadm_id,l.charttime,l.valuenum FROM mimiciv_hosp.labevents l
 JOIN project_d2v2.d2v2_map_lactate m ON l.itemid=m.itemid AND m.include_main
 WHERE l.valuenum IS NOT NULL AND l.valuenum>=0 AND EXISTS(SELECT 1 FROM base b WHERE b.hadm_id=l.hadm_id)
), lac AS (
 SELECT b.variant,b.stay_id,b.episode,max(l.valuenum) AS peak_lactate
 FROM base b LEFT JOIN lac_src l ON l.hadm_id=b.hadm_id
 AND l.charttime>=b.episode_start-interval '6 hours' AND l.charttime<=least(b.episode_start+interval '6 hours',b.t0)
 GROUP BY b.variant,b.stay_id,b.episode
), sep AS (
 SELECT b.variant,b.stay_id,b.episode,count(*) FILTER(WHERE s.sepsis3
 AND s.suspected_infection_time>=b.episode_start-interval '72 hours'
 AND s.suspected_infection_time<=least(b.episode_start+interval '24 hours',b.t0)) AS n_sepsis
 FROM base b LEFT JOIN mimiciv_derived.sepsis3 s ON s.stay_id=b.stay_id
 GROUP BY b.variant,b.stay_id,b.episode
)
SELECT b.*,i.inotrope,p.prbc_ml,l.peak_lactate,s.n_sepsis
FROM base b JOIN ino i USING(variant,stay_id,episode) JOIN prbc p USING(variant,stay_id,episode)
JOIN lac l USING(variant,stay_id,episode) JOIN sep s USING(variant,stay_id,episode)
ORDER BY variant,stay_id,t0
) TO STDOUT WITH CSV HEADER;
\o
\o :code_status_file
COPY (SELECT hadm_id,poe_id,ordertime,code_status FROM project_step4.step4_poe_code_status_orders ORDER BY hadm_id,ordertime,poe_id) TO STDOUT WITH CSV HEADER;
\o
\o :subject_last_record_file
COPY (SELECT subject_id,last_documented_dischtime FROM project_d4.d4_subject_last_record) TO STDOUT WITH CSV HEADER;
\o
ROLLBACK;
