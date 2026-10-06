"""Generate bounded, READ ONLY source query. Identifiers remain restricted."""
from pathlib import Path
import csv, json
from datetime import datetime

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
OUT=private_output('drug_scope');STAGE=OUT
source=release_path('mimic_selected_cohort')
with source.open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 1958 and len({r['stay_id'] for r in rows}) == 1958
target = STAGE / 'mimic_drug_windows.sql'
if target.exists():
    raise RuntimeError('refuse_query_overwrite')
values = []
for r in rows:
    key = int(r['stay_id'])
    q = datetime.fromisoformat(r['qualification_time']).isoformat(' ')
    h = datetime.fromisoformat(r['horizon']).isoformat(' ')
    values.append(f"({key},TIMESTAMP '{q}',TIMESTAMP '{h}')")
csvpath=(STAGE/'mimic_drug_windows.csv').resolve().as_posix()
query = f"""WITH fixed(stay_id,q,h) AS (VALUES {','.join(values)}) SELECT f.stay_id,v.starttime,v.endtime,v.norepinephrine,v.epinephrine,v.dopamine,v.phenylephrine,v.vasopressin FROM fixed f JOIN mimiciv_derived.vasoactive_agent v ON v.stay_id=f.stay_id AND v.starttime>=f.q AND v.starttime<=f.h WHERE v.endtime>v.starttime AND (v.norepinephrine>0 OR v.epinephrine>0 OR v.dopamine>0 OR v.phenylephrine>0 OR v.vasopressin>0) ORDER BY f.stay_id,v.starttime,v.endtime"""
target.write_text("\\set ON_ERROR_STOP on\nBEGIN TRANSACTION ISOLATION LEVEL REPEATABLE READ READ ONLY;\nSET LOCAL statement_timeout = '300s';\n\\copy (" + query + ") TO '" + csvpath + "' WITH (FORMAT csv, HEADER true)\nROLLBACK;\n", encoding='utf-8')
(OUT / 'mimic_extract_scope.json').write_text(json.dumps(dict(fixed_stays=1958,source_table='mimiciv_derived.vasoactive_agent',transaction='REPEATABLE READ READ ONLY, ROLLBACK',start_window='q <= start <= H',positive_rule='any of five recorded drug rates > 0; end > start',query_contains_identifiers_and_is_restricted=True), indent=2), encoding='utf-8')
print('Generated read-only query for 1958 fixed stays; no row data printed.')
