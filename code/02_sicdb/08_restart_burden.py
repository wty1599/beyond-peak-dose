"""S5 descriptive completion under the explicit overlap-UNKNOWN amendment.
No dose adjudication. All row-level data stay in restricted_work.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime
import csv, json, sqlite3, hashlib, math, time, sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
O=private_output('sicdb_burden');W=O/'private';W.mkdir(exist_ok=True)
start=time.monotonic();stamp=datetime.now().astimezone().isoformat()
src=release_path('sicdb_clinical_input');db=release_path('sicdb_sqlite')
def sha(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def csvout(path,rows,fields=None):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0])); w.writeheader(); w.writerows(rows)
def round_nee(x): return math.floor(x*10000+.5)/10000
before={str(p):sha(p) for p in (src,db)}
csvout(O/'S5_resumed_source_before.csv',[dict(path=p,sha256=h) for p,h in before.items()])
with src.open(encoding='utf-8-sig',newline='') as f: allrows=list(csv.DictReader(f))
rows=[r for r in allrows if r['restart_binary'] in ('1','1.0')]
assert len(allrows)==2044 and len(rows)==684
c=sqlite3.connect(db.resolve().as_uri()+'?mode=ro',uri=True); c.execute('PRAGMA query_only=ON')
assert c.execute('PRAGMA query_only').fetchone()[0]==1
F={1562:1.,1502:1.,1618:.01,1593:.1,1550:2.5}
record=[]; ledger=[]; overlap_log=[]; failures=[]; raw_count=0
for r in rows:
    try:
        cid=int(r['CaseID']); first=float(r['valid_restart_time']); end=float(r['ascertained_followup_end']); t0=float(r['t0'])
        assert first>=t0+14400 and first<end<=t0+72*3600
        weight=c.execute('SELECT WeightOnAdmission FROM cases WHERE CaseID=?',(cid,)).fetchone()[0]/1000
        assert math.isfinite(weight) and weight>0
        meds=c.execute('SELECT id,DrugID,Offset,OffsetDrugEnd,AmountPerMinute FROM medication WHERE CaseID=? AND DrugID IN (1562,1502,1618,1593,1550) AND AmountPerMinute>0 AND OffsetDrugEnd>Offset AND Offset<? AND OffsetDrugEnd>? ORDER BY Offset,id',(cid,end,first)).fetchall()
        raw_count+=len(meds); edges=defaultdict(list)
        for mid,drug,a,b,rate in meds:
            assert b>a and math.isfinite(rate)
            dose=rate*F[drug]
            if drug!=1550: dose=dose*1e6/weight
            assert math.isfinite(dose) and dose>0
            edges[max(float(a),first)].append((1,mid,drug,dose)); edges[float(b)].append((-1,mid,drug,dose))
        times=sorted(edges); active={k:{} for k in F}; segments=[]; ambiguous=False
        for i,a in enumerate(times[:-1]):
            for direction,mid,drug,dose in edges[a]:
                if direction==1: active[drug][mid]=dose
                else: active[drug].pop(mid,None)
            if a>=end: break
            used=[k for k in F if active[k]]
            if not used: continue
            b=times[i+1]; cb=min(b,end)
            if cb<=a: continue
            amb=any(len(active[k])>1 for k in used)
            if amb:
                ambiguous=True; nee=None
                # A positive single non-overlapped drug proves positivity without
                # resolving the other doses. Alternatively every possible active
                # source record is itself positive at the original precision.
                independently_positive=(any(len(active[k])==1 and round_nee(next(iter(active[k].values())))>0 for k in used)
                                        or all(round_nee(d)>0 for k in used for d in active[k].values()))
                pos=1 if independently_positive else None
                overlap_log.append(dict(CaseID=cid,start=a,end=cb,positivity='POSITIVE_INDEPENDENTLY_KNOWN' if pos==1 else 'UNKNOWN'))
            else:
                nee=round_nee(sum(next(iter(active[k].values())) for k in used)); pos=int(nee>0)
            seg=dict(CaseID=cid,raw_start=a,raw_end=b,clipped_end=cb,nee=nee,positive_state=pos,same_drug_overlap=amb)
            ledger.append(seg); segments.append(seg)
        known=[s for s in segments if s['positive_state']==1]
        uncertain=[s for s in segments if s['positive_state'] is None]
        assert known and known[0]['raw_start']==first
        hours=None if uncertain else sum((s['clipped_end']-s['raw_start'])/3600 for s in known)
        if hours is not None: assert 0<hours<=(end-first)/3600+1e-10
        def category(last):
            return 'POSITIVE_RECORD_CROSSES_OBSERVATION_END' if last>end else 'TERMINAL_TIE_COMPLETENESS_UNRESOLVED' if last==end else 'CESSATION_OBSERVED_BEFORE_WINDOW_END'
        latest_known=max(s['raw_end'] for s in known)
        cls=category(latest_known)
        if uncertain:
            possible_latest=max(latest_known,max(s['raw_end'] for s in uncertain))
            if category(possible_latest)!=cls: cls='UNKNOWN_DUE_TO_AMBIGUOUS_POSITIVITY'
        peak=None if ambiguous else max(s['nee'] for s in known)
        record.append(dict(CaseID=cid,PatientID=int(r['PatientID']),first_restart=first,followup_end=end,
            end_reason=r['followup_end_reason'],completeness=cls,positive_hours=hours,peak_nee=peak,
            initial_episode_peak_nee=float(r['predictor_episode_peak_nee']),observable_hours=(end-first)/3600,
            hours_to_restart=(first-t0)/3600,same_drug_overlap=ambiguous,positive_state_unknown_segments=len(uncertain),
            peak_status='UNKNOWN_SAME_DRUG_OVERLAP' if ambiguous else 'KNOWN',
            duration_status='UNKNOWN_POSITIVITY' if uncertain else 'KNOWN'))
    except Exception as exc:
        failures.append(dict(CaseID=int(r['CaseID']),reason=str(exc)))
        csvout(W/'S5_resumed_failures.csv',failures)
        csvout(O/'S5_resumed_failures_aggregate.csv',[dict(failure_count=len(failures),reason=str(exc))])
        c.close(); raise
c.close()
csvout(W/'S5_positive_nee_segments.csv',ledger)
csvout(W/'S5_burden_event_records.csv',record)
csvout(W/'S5_same_drug_overlap_segments.csv',overlap_log,['CaseID','start','end','positivity'])
summary=[]
for group in ['ALL']+sorted({r['completeness'] for r in record}):
    sub=[r for r in record if group=='ALL' or r['completeness']==group]
    for measure in ['positive_hours','peak_nee','initial_episode_peak_nee','observable_hours']:
        x=np.array([r[measure] for r in sub if r[measure] is not None],dtype=float)
        assert np.isfinite(x).all()
        q=np.quantile(x,[.25,.5,.75],method='linear') if len(x) else [None]*3
        summary.append(dict(group=group,measure=measure,n=len(sub),valid_n=len(x),unknown_n=len(sub)-len(x),
            q1=q[0],median=q[1],q3=q[2],role='observation metadata' if measure=='observable_hours' else 'descriptive treatment burden or initial exposure'))
csvout(O/'S5_restart_burden_summary.csv',summary)
counts=defaultdict(int)
for r in record: counts[(r['end_reason'],r['completeness'])]+=1
csvout(O/'S5_observation_end_and_completeness.csv',[dict(end_reason=k[0],completeness=k[1],cases=n) for k,n in sorted(counts.items())])
availability=[dict(measure=m,total_n=684,known_n=sum(r[m] is not None for r in record),unknown_n=sum(r[m] is None for r in record)) for m in ['positive_hours','peak_nee','initial_episode_peak_nee']]
csvout(O/'S5_measure_availability.csv',availability)
csvout(O/'S5_resumed_failures_aggregate.csv',[],['failure_count','reason'])
after={str(p):sha(p) for p in (src,db)}; assert before==after
csvout(O/'S5_resumed_source_hash_check.csv',[dict(path=p,before=h,after=after[p],unchanged=True) for p,h in before.items()])
qc=dict(status='COMPLETE_AWAITING_INDEPENDENT_REVIEW',cases=len(record),patients=len({r['PatientID'] for r in record}),
    same_drug_overlap_cases=sum(r['same_drug_overlap'] for r in record),overlap_segments=len(overlap_log),
    positivity_unknown_segments=sum(r['positive_state_unknown_segments'] for r in record),
    duration_unknown=sum(r['positive_hours'] is None for r in record),peak_nee_unknown=sum(r['peak_nee'] is None for r in record),
    completeness_unknown=sum(r['completeness']=='UNKNOWN_DUE_TO_AMBIGUOUS_POSITIVITY' for r in record),
    source_medication_records=raw_count,source_unchanged=True,sqlite_read_only=True,author_amendment_applied=True,
    seed='not applicable; deterministic description',python=sys.version,numpy=np.__version__,
    elapsed_seconds=time.monotonic()-start,started_at=stamp,ended_at=datetime.now().astimezone().isoformat(),failures=0)
(O/'S5_completed_QC.json').write_text(json.dumps(qc,indent=2),encoding='utf-8')
print(json.dumps(qc,indent=2));print(json.dumps([r for r in summary if r['group']=='ALL'],indent=2))
