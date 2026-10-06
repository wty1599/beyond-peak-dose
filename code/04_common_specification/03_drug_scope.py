"""Fixed-cohort five/four-drug endpoint bridge; no cohort reconstruction."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timedelta
import csv, gzip, hashlib, json, math, sqlite3, sys, time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
OUT=private_output('drug_scope');STAGE=OUT
DRUGS=('norepinephrine','epinephrine','dopamine','phenylephrine','vasopressin')
SIDS={1562:'norepinephrine',1502:'epinephrine',1618:'dopamine',1593:'phenylephrine',1550:'vasopressin'}

def read(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

def write(path,rows,fields=None):
    if path.exists(): raise RuntimeError('refuse_output_overwrite:'+path.name)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)

def num(x):
    try: y=float(x)
    except (TypeError,ValueError): return None
    return y if math.isfinite(y) else None

def dt(x): return datetime.fromisoformat(x) if x and x not in ('NaT','nan','None') else None
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def valid_start(a,q,h,e,d=None,coverage=None,day=None):
    return (q<=a<=h and a<e and (d is None or a<d) and
            (coverage is None or a<coverage) and (day is None or d is not None or a<day))

def boundary_tests():
    assert valid_start(4,4,72,100)
    assert valid_start(72,4,72,100)
    assert not valid_start(10,4,72,10)
    assert not valid_start(10,4,72,100,d=10)
    assert not valid_start(10,4,72,100,coverage=10)
    # First Vaso-only does not remove a later four-drug start.
    rows=[(4,{'vasopressin'}),(8,{'norepinephrine'})]
    assert min(t for t,ds in rows if ds-{'vasopressin'})==8
    assert not valid_start(10,4,72,100,day=9)
    return 7

def mimic():
    cohort=release_path('mimic_selected_cohort')
    drugfile=STAGE/'mimic_drug_windows.csv'
    rows=read(cohort); raw=read(drugfile)
    assert len(rows)==1958
    group=defaultdict(list)
    for z in raw:
        ds={d for d in DRUGS if (num(z[d]) or 0)>0}
        assert ds and dt(z['endtime'])>dt(z['starttime'])
        group[int(z['stay_id'])].append((dt(z['starttime']),ds))
    outputs=[]
    for r in rows:
        q,h,e=dt(r['qualification_time']),dt(r['horizon']),dt(r['care_end'])
        d,day=dt(r['deathtime']),dt(r['dod'])
        intervals=group[int(r['stay_id'])]
        accepted=[(a,ds) for a,ds in intervals if valid_start(a,q,h,e,d=d,day=day)]
        five=bool(accepted)
        original=int(float(r['restart_72h']))
        if int(five)!=original: raise RuntimeError('MIMIC_five_drug_label_parity_failed; row remains restricted')
        first=min((a for a,_ in accepted),default=None)
        first_set=set().union(*(ds for a,ds in accepted if a==first)) if first else set()
        four_starts=[a for a,ds in accepted if ds-{'vasopressin'}]
        four=1 if four_starts else 0
        # Reuse the original date-level ordering policy for a new first four-drug start.
        # All main-cohort ascertained ends were checked before extraction; still keep this branch explicit.
        ambiguous=[a for a,ds in intervals if ds-{'vasopressin'} and q<=a<=h and a<e and d is None and day is not None and day<=a<day+timedelta(days=1)]
        if not four and ambiguous: four=None
        original_time=dt(r['effective_next_positive'])
        if five and first!=original_time: raise RuntimeError('MIMIC_first_restart_time_parity_failed')
        outputs.append(dict(database='MIMIC-IV',event_key=r['stay_id'],patient_key=r['subject_id'],five_label=original,four_label=four,
            first_vaso_only=int(first_set=={'vasopressin'}),first_vaso_then_four=int(first_set=={'vasopressin'} and four==1),
            window_vaso_only=int(five and four==0),first_five_time=first.isoformat() if first else '',first_four_time=min(four_starts).isoformat() if four_starts else ''))
    return outputs,dict(source_cohort_sha256=sha(cohort),drug_extract_sha256=sha(drugfile),extracted_positive_intervals=len(raw))

def sicdb():
    cohort=release_path('sicdb_clinical_input')
    db=release_path('sicdb_sqlite')
    allrows=read(cohort); assert len(allrows)==2044
    rows=[r for r in allrows if num(r['restart_binary']) is not None]
    assert sum(int(float(r['restart_binary'])) for r in rows)==684
    con=sqlite3.connect(db.resolve().as_uri()+'?mode=ro',uri=True);con.execute('PRAGMA query_only=ON')
    outputs=[];count=0
    try:
        for r in rows:
            q,h,e=num(r['q']),num(r['horizon']),num(r['care_end'])
            d,cov=num(r['death']),num(r['coverage_end'])
            assert q is not None and h is not None and e is not None and cov is not None
            raw=con.execute('SELECT DrugID,Offset,OffsetDrugEnd,AmountPerMinute FROM medication WHERE CaseID=? AND DrugID IN (1562,1502,1618,1593,1550) AND AmountPerMinute>0 AND OffsetDrugEnd>Offset AND Offset>=? AND Offset<=? ORDER BY Offset,DrugID',(int(r['CaseID']),q,h)).fetchall()
            count+=len(raw)
            accepted=[(float(a),SIDS[did]) for did,a,b,v in raw if valid_start(float(a),q,h,e,d=d,coverage=cov)]
            original=int(float(r['restart_binary']))
            if int(bool(accepted))!=original: raise RuntimeError('SICdb_five_drug_label_parity_failed; row remains restricted')
            first=min((a for a,_ in accepted),default=None)
            first_set={drug for a,drug in accepted if a==first}
            four_starts=[a for a,drug in accepted if drug!='vasopressin']
            four=1 if four_starts else None
            if four is None:
                firstD=d is not None and q<d<=h and d<=e and d<=cov
                normal=r['closure_class']=='NORMAL_CARE_END' and r['coverage_type']=='FULL_MEDICATION_TABLE_AND_NORMAL_CARE_END'
                observed_h=cov>=h and e>=h
                if firstD or normal or observed_h: four=0
            if original==0 and four is None:
                # Original zero establishes a definite observed no-drug endpoint for this same immutable window.
                four=0
            original_time=num(r['valid_restart_time'])
            if original==1 and (first is None or first!=original_time): raise RuntimeError('SICdb_first_restart_time_parity_failed')
            outputs.append(dict(database='SICdb',event_key=r['CaseID'],patient_key=r['PatientID'],five_label=original,four_label=four,
                first_vaso_only=int(first_set=={'vasopressin'}),first_vaso_then_four=int(first_set=={'vasopressin'} and four==1),
                window_vaso_only=int(original==1 and four==0),first_five_time=first if first is not None else '',first_four_time=min(four_starts) if four_starts else ''))
    finally: con.close()
    return outputs,dict(source_cohort_sha256=sha(cohort),sqlite_sha256=sha(db),extracted_positive_intervals=count,source_total_cases=len(allrows),excluded_original_restart_unknown=len(allrows)-len(rows))

def main():
    started=time.time(); db=sys.argv[1]
    assert db in ('mimic','sicdb')
    tests=boundary_tests()
    rows,audit=(mimic if db=='mimic' else sicdb)()
    assert all(r['four_label']!=1 or r['five_label']==1 for r in rows)
    known=[r for r in rows if r['four_label'] is not None]
    n=len(rows);nf=sum(r['five_label'] for r in rows);n4=sum(r['four_label'] for r in known)
    summary=dict(database=rows[0]['database'],fixed_events=n,patients=len({r['patient_key'] for r in rows}),five_restarts=nf,five_evaluable=n,five_percent=100*nf/n,
      four_restarts=n4,four_evaluable=len(known),four_unknown=n-len(known),four_percent=100*n4/len(known),
      first_restart_vaso_only=sum(r['first_vaso_only'] for r in rows),first_vaso_then_four=sum(r['first_vaso_then_four'] for r in rows),
      valid_window_vaso_only=sum(r['window_vaso_only'] for r in rows),
      paired_restart_count_difference=n4-sum(r['five_label'] for r in known),
      paired_percentage_point_difference=100*(n4-sum(r['five_label'] for r in known))/len(known),
      difference_denominator='same four-drug-evaluable fixed events',new_cohort_reconstruction=False)
    write(STAGE/(db+'_bridge_event_ledger.csv'),rows)
    write(OUT/(db+'_fixed_cohort_outcome_bridge.csv'),[summary])
    audit.update(boundary_tests=tests,five_label_mismatches=0,first_restart_time_mismatches=0,full_window_scanned=True,no_NEE_reestimated=True,no_main_cohort_modified=True,no_public_identifiers=True,elapsed_seconds=time.time()-started)
    (OUT/(db+'_bridge_audit.json')).write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False))

if __name__=='__main__': main()
