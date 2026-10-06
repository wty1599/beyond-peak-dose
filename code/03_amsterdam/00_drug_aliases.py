"""Targeted names and candidate-unit scan only. No cohort or outcome calculation."""
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime, timezone, timedelta
import csv, gzip, hashlib, json, math, re, time

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
OUT=private_output('amsterdam_alias');STAGE=OUT
SOURCE=release_path('amsterdam_raw')/'drugitems.csv'
EXPECTED='3bf8bdbabf3e7b67d7ad6d4b9389f2fd689d70d28c6052cb8013ad4a3cddcaa9'
EXPECTED_ROWS=4907269
TERMS=('vasopressin','vasopressine','argipressin','argipressine','pitressin',
       'empressin','empressine','vasostrict','antidiuret')
ANALOG=('terlipressin','terlipressine','glypressin')
target=OUT/'Amsterdam_drug_alias_check.json'
if target.exists():raise RuntimeError('refuse_overwrite')
assert SOURCE.stat().st_size==818558125
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
assert sha(SOURCE)==EXPECTED
started=time.monotonic()
counts=Counter();per_drug=Counter();admissions=defaultdict(set);units=Counter();candidate_rows=[]
csv.field_size_limit(16*1024*1024)
with SOURCE.open(encoding='cp1252',newline='') as f:
    reader=csv.DictReader(f)
    columns=reader.fieldnames
    assert all(k in columns for k in ['item','solutionitem','admissionid','itemid','start','stop','rate','administered'])
    n=0
    for row in reader:
        n+=1
        if n>EXPECTED_ROWS or None in row:raise RuntimeError('source_structure_changed')
        names=' '.join([row['item'],row['solutionitem']]).casefold()
        is_analog=any(t in names for t in ANALOG)
        is_vaso=any(t in names for t in TERMS) or re.search(r'\badh\b',names) is not None
        if not (is_analog or is_vaso):continue
        group='terlipressin_analog_not_in_NEE' if is_analog else 'vasopressin_name_candidate'
        counts[group]+=1
        key=(group,row['itemid'],row['item'],row['solutionitem'])
        per_drug[key]+=1;admissions[group].add(row['admissionid'])
        units[(group,row['doseunit'],row['doserateunit'],row['administeredunit'],row['rateunit'])]+=1
        candidate_rows.append({**row,'candidate_group':group})
assert n==EXPECTED_ROWS and sha(SOURCE)==EXPECTED
private=STAGE/'Amsterdam_alias_candidates.csv.gz'
with gzip.open(private,'wt',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=columns+['candidate_group']);w.writeheader();w.writerows(candidate_rows)
# Do not release complementary name/solution subcounts alongside the total.
# Only frequently represented, public dictionary identities are retained.
named=sorted({(k[0],k[1],k[2]) for k,value in per_drug.items() if value>=5})
result=dict(execution_time=datetime.now(timezone(timedelta(hours=8))).isoformat(),
            source=str(SOURCE),source_sha256=EXPECTED,source_rows=n,
            searched_fields=['item','solutionitem'],vasopressin_terms=list(TERMS)+['word_ADH'],
            analog_terms=list(ANALOG),
            vaso_candidate_rows=counts['vasopressin_name_candidate'],
            vaso_candidate_admissions=len(admissions['vasopressin_name_candidate']),
            analog_order_rows=counts['terlipressin_analog_not_in_NEE'],
            analog_distinct_admissions=len(admissions['terlipressin_analog_not_in_NEE']),
            mapped_names_without_subcounts=[dict(group=k[0],itemid=k[1],item=k[2]) for k in named],
            solution_subcounts_not_released=True,
            source_unchanged=True,no_outcomes_accessed=True,no_cohort_modified=True,
            candidate_rows_restricted=True,elapsed_seconds=time.monotonic()-started,
            conclusion='Names negative does not prove clinical absence; terlipressin is excluded from NEE.')
target.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:result[k] for k in ('source_rows','vaso_candidate_rows','analog_order_rows',
                                       'analog_distinct_admissions','source_unchanged','no_outcomes_accessed')},ensure_ascii=False))
