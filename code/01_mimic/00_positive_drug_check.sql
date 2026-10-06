\set ON_ERROR_STOP on
\pset pager off
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout='600s';
\o :positive_drug_file
COPY (
 SELECT v.stay_id,v.starttime,v.endtime,n.stay_id IS NULL AS missing_nee_interval
 FROM mimiciv_derived.vasoactive_agent v
 LEFT JOIN mimiciv_derived.norepinephrine_equivalent_dose n
 ON n.stay_id=v.stay_id AND n.starttime=v.starttime AND n.endtime=v.endtime
 WHERE (v.norepinephrine>0 OR v.epinephrine>0 OR v.dopamine>0 OR v.phenylephrine>0 OR v.vasopressin>0)
 AND coalesce(n.norepinephrine_equivalent_dose,0)<=0
 ORDER BY v.stay_id,v.starttime,v.endtime
) TO STDOUT WITH CSV HEADER;
\o
ROLLBACK;
