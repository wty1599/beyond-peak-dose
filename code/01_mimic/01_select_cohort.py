"""Build the final cohort from candidate and source-evidence exports."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timedelta
import csv
import importlib.util
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
spec=importlib.util.spec_from_file_location('definitions',Path(__file__).with_name('01_outcome_definitions.py'))
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
dt, normalise, terminal_conflicts, qualification, evaluate=d.dt,d.normalise,d.terminal_conflicts,d.qualification,d.evaluate
Q=d.Q
INPUT=release_path('mimic_candidate_exports')
OUTPUT=private_output('mimic_cohort')
if (OUTPUT/'selected_cohort.csv').exists():raise FileExistsError('Cohort output already exists.')
def read(name):
    with (INPUT/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUTPUT/name).open('x',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
candidates=[normalise(r) for r in read('candidates.csv')]
zero_drugs=defaultdict(list)
for z in read('positive_drug_rounded_zero_nee.csv'):
    zero_drugs[int(z['stay_id'])].append((dt(z['starttime']),dt(z['endtime'])))
for r in candidates:
    intervals=zero_drugs[r['stay_id']]
    r['rounded_zero_drug_in_qualification']=any(a<r['t0']+Q and b>r['t0'] for a,b in intervals)
    later=[a for a,b in intervals if a>=r['t0'] and a<r['outtime']]
    r['effective_next_positive']=min(later+([r['next_positive']] if r['next_positive'] else [])) if later or r['next_positive'] else None
    r['terminal_source_conflicts']=';'.join(terminal_conflicts(r))
orders=defaultdict(list)
for r in read('code_status.csv'):orders[int(r['hadm_id'])].append((dt(r['ordertime']),r['poe_id'],r['code_status']))
for a in orders.values():a.sort()
last={int(r['subject_id']):dt(r['last_documented_dischtime']) for r in read('subject_last_record.csv')}
def lim(r):
    oo=orders[r['hadm_id']]; prior=[v for v in oo if v[0]<r['t0']]
    if prior and prior[-1][2]=='LIMITATION_PRESENT':return 'LIMITATION_PRE'
    if any(x[0]>=r['t0'] and x[2]=='LIMITATION_PRESENT' for x in oo) or 'hospice' in r['discharge_location'].lower():return 'LIMITATION_POST'
    return 'LIMITATION_NONE' if oo else 'LIMITATION_UNKNOWN'
flows=[]; groups=defaultdict(list)
for r in candidates:
    r['clinical_pass']=not r['inotrope'] and r['prbc_ml']<1400 and r['peak_lactate'] is not None and r['peak_lactate']>2 and r['n_sepsis']>0
    r['limitation_state']=lim(r)
    r['time_status'],r['time_reason']=qualification(r)
    groups[(r['variant'],r['stay_id'])].append(r)
selected=[];selection_status={};unknown_stays=[];pre_stays=[]
for (v,s),rr in groups.items():
    if v!='new':continue
    rr=sorted([r for r in rr if r['clinical_pass']],key=lambda r:r['t0'])
    selection_status[s]='NO_TIME_QUALIFIED_CLINICAL_CANDIDATE'
    for r in rr:
        if r['time_status']=='FAIL':continue
        if r['time_status']=='UNKNOWN':
            selection_status[s]='EARLIER_TIME_UNASSESSABLE';unknown_stays.append(r);break
        if r['limitation_state']=='LIMITATION_PRE':
            selection_status[s]='FIRST_QUALIFIED_PRE_EXCLUDED';pre_stays.append(r);break
        selection_status[s]='SELECTED';r.update(evaluate(r,last.get(r['subject_id'])));selected.append(r);break

write('selected_cohort.csv',selected)
write('candidate_classifications.csv',candidates)
write('first_identity_unknown.csv',unknown_stays)
write('first_pre_exclusions.csv',pre_stays)
