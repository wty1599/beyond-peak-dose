"""Prepare a bounded, read-only extract for already-adjudicated restarts."""
from pathlib import Path
from datetime import datetime
import csv, hashlib, json, sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output

root=private_output('mimic_burden')
w=root/'restricted_work';w.mkdir(exist_ok=True)
p=root.parents[1]
source=release_path('mimic_selected_cohort')
with source.open(encoding='utf-8-sig',newline='') as f:
    rows=[r for r in csv.DictReader(f) if r['restart_72h']=='1']
assert len(rows)==680
before=hashlib.file_digest(source.open('rb'),'sha256').hexdigest()
values=[]
for r in rows:
    start=r['cif_time'];end=r['ascertained_followup_end']
    assert r['cif_status']=='RESTART' and start and end
    assert datetime.fromisoformat(start)<datetime.fromisoformat(end)
    assert start==r['effective_next_positive']
    values.append('('+r['stay_id']+",'"+start+"'::timestamp,'"+end+"'::timestamp)")
path=(w/'nee_after_restart.csv').resolve().as_posix()
sql="""\\set ON_ERROR_STOP on
\\pset pager off
BEGIN TRANSACTION ISOLATION LEVEL REPEATABLE READ READ ONLY;
SET LOCAL statement_timeout='120s';
SELECT current_setting('transaction_read_only') AS read_only;
\\o OUTPUT
COPY (
 WITH requested(stay_id,first_restart,followup_end) AS (VALUES VALUES_LIST)
 SELECT r.stay_id,n.starttime AS raw_starttime,n.endtime AS raw_endtime,
 n.norepinephrine_equivalent_dose::double precision AS nee
 FROM requested r JOIN mimiciv_derived.norepinephrine_equivalent_dose n USING(stay_id)
 WHERE n.norepinephrine_equivalent_dose>0
 AND n.endtime>r.first_restart AND n.starttime<r.followup_end
 ORDER BY r.stay_id,n.starttime,n.endtime
) TO STDOUT WITH CSV HEADER;
\\o
ROLLBACK;
""".replace('OUTPUT',path).replace('VALUES_LIST',',\n'.join(values))
(w/'extract_readonly.sql').write_text(sql,encoding='utf-8')
(root/'pre_extraction_record.json').write_text(json.dumps(dict(status='BEFORE_EXTRACTION',time=datetime.now().astimezone().isoformat(),n=680,source_sha256=before,
 only_existing_restarts=True,outcomes_changed=False,uncertain_followup_end=0,
 window='[first adjudicated restart, existing ascertained_followup_end)',
 duration='Union of positive total-NEE intervals clipped to window; convert each disjoint segment to hours then sum',
 peak='Maximum existing total NEE among positive intersections',
 complete_definition='Last positive record ends strictly before followup_end; describes within-window cessation only',
 truncation='Record crosses followup_end: truncated; record ends exactly at followup_end: terminal-tie unresolved, not complete',
 summaries='Median and Q1/Q3 with valid n, overall and by completeness; followup-window length only as observability metadata',
 no_thresholds_no_models_no_tests=True),indent=2),encoding='utf-8')
print('Prepared read-only bounded interval extraction for 680 existing restarts; no row identifiers printed.')
