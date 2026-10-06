"""Read approved inputs directly; calculate requested cells, not upstream labels."""
from pathlib import Path
from collections import Counter
from datetime import datetime
import csv, gzip, json, math, re, sqlite3

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
STAGE=private_output('amsterdam_descriptive');OUT=STAGE;ROOT=Path('.')
HOUR = 3600000

def read(path, encoding='utf-8-sig'):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding=encoding, newline='') as handle:
        return list(csv.DictReader(handle))

def write_private(name, rows):
    with (STAGE/name).open('w', encoding='utf-8', newline='') as handle:
        w = csv.DictWriter(handle,fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def num(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError,ValueError): return None

paths = []
def load(relative, encoding='utf-8-sig'):
    path = ROOT/relative
    paths.append(str(path))
    return read(path,encoding)

cohort=read(release_path('amsterdam_cohort_input'))
admissions_path = release_path('amsterdam_admissions')
paths.append(str(admissions_path))
admissions = {r['admissionid']:r for r in read(admissions_path,'cp1252')}
if len(cohort)!=779 or sum(r['h72_restart'] in ('YES','NO') for r in cohort)!=766:
    raise RuntimeError('STOP: Amsterdam denominator differs from 779/766')

def band(row):
    lo,hi = num(row['four_drug_nee_peak_lower']),num(row['four_drug_nee_peak_upper'])
    if hi is not None and hi<.20: return '<0.20'
    if lo is not None and hi is not None and lo>=.20 and hi<.50: return '0.20-<0.50'
    if lo is not None and lo>=.50: return '>=0.50'
    return 'unknown'

counts = []
a1_rows = []
for r in cohort:
    a = admissions[r['admissionid']]
    landmark = int(r['t0_ms'])+24*HOUR
    care_end = num(a['dischargedat'])
    # Three-valued conjunction. Definite restart/end/death excludes; uncertain
    # survival, same-unit observation or restart-free state remains unknown.
    if r['h24_restart']=='YES' or r['h24_death']=='YES' or (care_end is not None and care_end<=landmark):
        condition='excluded'
    elif care_end is None or r['h24_restart']!='NO' or r['h24_death']!='NO':
        condition='unknown'
    else:
        condition='eligible'
    a1_rows.append(dict(patientid=r['patientid'],condition=condition,
                        restart72=r['h72_restart'],death24=r['h24_death']))
write_private('a1_rows.csv',a1_rows)
eligible=[r for r in a1_rows if r['condition']=='eligible']
unknown=[r for r in a1_rows if r['condition']=='unknown']
known=[r for r in eligible if r['restart72'] in ('YES','NO')]
a1=dict(condition_n=len(eligible),restart_n=sum(r['restart72']=='YES' for r in known),
        outcome_N=len(known),outcome_unknown=len(eligible)-len(known),state_unknown=len(unknown),
        state_unknown_date=sum(r['death24']=='UNKNOWN_DATE_PRECISION' for r in unknown),
        excluded=sum(r['condition']=='excluded' for r in a1_rows))
# Requested extreme label range: include indeterminate conditional-status cases
# as well as any indeterminate future labels. This is a sensitivity envelope,
# not an estimated eligible population or confidence interval.
potential=len(eligible)+len(unknown)
a1['unknown_negative_bound_percent']=100*a1['restart_n']/potential
a1['unknown_positive_bound_percent']=100*(a1['restart_n']+len(unknown)+a1['outcome_unknown'])/potential
a1['bound_N']=potential

a2=[]
det=[r for r in cohort if r['h72_restart'] in ('YES','NO')]
for variable,field in (('Ventilation process','recorded_ventilation_process_at_t0'),
                       ('CVVH process','recorded_cvvh_process_at_t0')):
    for level,state in (('Recorded','RECORDED_ACTIVE_PROCESS'),('Not recorded','NO_MAPPED_PROCESS_RECORDED')):
        rows=[r for r in det if r[field]==state]
        a2.append(dict(variable=variable,level=level,n=sum(r['h72_restart']=='YES' for r in rows),N=len(rows)))
    unknown_count=sum(r[field] not in ('RECORDED_ACTIVE_PROCESS','NO_MAPPED_PROCESS_RECORDED') for r in det)
    a2.append(dict(variable=variable,level='unknown',n=unknown_count,N=len(det)))
for level in ('<0.20','0.20-<0.50','>=0.50','unknown'):
    rows=[r for r in det if band(r)==level]
    a2.append(dict(variable='Peak four-drug NEE band',level=level,
                   n=(len(rows) if level=='unknown' else sum(r['h72_restart']=='YES' for r in rows)),
                   N=(len(det) if level=='unknown' else len(rows))))

def dist(values): return dict(Counter(values))
a4=dict(age=dist(admissions[r['admissionid']]['agegroup'] or 'unknown' for r in cohort),
        sex=dist(admissions[r['admissionid']]['gender'] or 'unknown' for r in cohort),
        peak=dist(band(r) for r in cohort),
        ventilation=dist(r['recorded_ventilation_process_at_t0'] for r in cohort),
        cvvh=dist(r['recorded_cvvh_process_at_t0'] for r in cohort))
write_private('a4_duration.csv',[dict(hours=r['positive_hours_segmentwise']) for r in cohort])

alias=read(release_path('amsterdam_alias_input'))
cohort_by_admission={r['admissionid']:r for r in cohort}
a5_ids=set(); a5_time_unknown=0
for r in alias:
    if r['admissionid'] not in cohort_by_admission: continue
    if not re.search(r'terlipress|glypressin',r.get('item','')+' '+r.get('solutionitem',''),re.I): continue
    t0=int(cohort_by_admission[r['admissionid']]['t0_ms'])
    start,stop=num(r['start']),num(r['stop'])
    if start is None: a5_time_unknown+=1; continue
    low,high=t0-24*HOUR,t0+72*HOUR
    if low<=start<=high or (stop is not None and start<low<stop): a5_ids.add(r['admissionid'])
a5=dict(events=len(a5_ids),N=779,time_unknown_orders=a5_time_unknown)


summary=dict(A1=a1,A2=a2,A4=a4,A5=a5)
(STAGE/'descriptive_cells.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
