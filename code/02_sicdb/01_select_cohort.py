"""Approved chronology-only SICdb rebuild. New outputs, no models/source writes."""
from pathlib import Path
from collections import defaultdict, Counter
import argparse, hashlib, json, math, sqlite3, sys, time
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
RAW=release_path('sicdb_raw')
WORK=private_output('sicdb_cohort');AGG=WORK
if (WORK/'selected_cohort.csv.gz').exists():raise FileExistsError('Cohort output already exists.')
FIVE={1562,1502,1618,1593,1550}
FACTOR={1562:1.,1502:1.,1618:.01,1593:.1,1550:2.5}
HOSP_ALIVE={2026,3129,3131,3132,3133,3134}; HOSP_DEAD={2028,3130}
NA=['','NULL','null','\\N','NA','NaN']
Q=14400; H=72*3600
QC=Counter()
def log(x): print(time.strftime('%H:%M:%S'),x,flush=True)
def finite(x): return x is not None and pd.notna(x) and np.isfinite(x)
def jwrite(p,x): p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)),encoding='utf-8')
def merge_ranges(ranges):
    out=[]
    for a,b in sorted(ranges):
        if not out or a>out[-1][1]:out.append([a,b])
        else:out[-1][1]=max(out[-1][1],b)
    return out

cases=pd.read_csv(RAW/'cases.csv.gz',na_values=NA).set_index('CaseID',drop=False)
valid=cases[(cases.AgeOnAdmission>=18)&(cases.TimeOfStay>cases.ICUOffset)&(cases.WeightOnAdmission>0)]
ids=set(valid.index)
# Optional evidence for end-unknown cases; no threshold/frequency criterion.
# Include every non-alive-discharge case; hospital type alone may lack a day.
unclear=set(valid.index[~valid.DischargeState.eq(2202)])
extra=defaultdict(list); parts=[]; total=0
columns=['id','CaseID','DrugID','Offset','OffsetDrugEnd','Amount','AmountPerMinute','IsSingleDose','GivenState']
for ch in pd.read_csv(RAW/'medication.csv.gz',na_values=NA,chunksize=250000,usecols=columns):
    total+=len(ch)
    keep=ch[ch.CaseID.isin(ids)&ch.DrugID.isin(FIVE|{1559,1560,2046})]
    if len(keep):parts.append(keep)
    oth=ch[ch.CaseID.isin(unclear)&~ch.DrugID.isin(FIVE)&(ch.AmountPerMinute>0)&(ch.OffsetDrugEnd>ch.Offset)]
    for r in oth.itertuples(index=False):extra[int(r.CaseID)].append((float(r.Offset),float(r.OffsetDrugEnd)))
med=pd.concat(parts,ignore_index=True)
extra={k:merge_ranges(v) for k,v in extra.items()}
QC['raw_medication_rows_read']=total
assert total==5141127,'Raw table incomplete or changed; stop'
log('Full medication source read completed')
con=sqlite3.connect(release_path('sicdb_sqlite').resolve().as_uri()+'?mode=ro',uri=True)
con.execute('PRAGMA query_only=ON')
ecmo=set(x[0] for x in con.execute('SELECT DISTINCT CaseID FROM signals WHERE DataID IN (2023,2024,2025) AND Val>0'))
labs=pd.read_sql_query('SELECT CaseID,Offset,LaboratoryValue FROM labs WHERE LaboratoryID IN (454,465,657)',con)
labs=labs[np.isfinite(labs.LaboratoryValue)&(labs.LaboratoryValue>=0)]
labgroups={int(k):v for k,v in labs.groupby('CaseID',sort=False)}
con.close()
mapping=pd.read_csv(release_path('shock_icd_mapping'))
mapping=mapping[(mapping.icd_version==10)&mapping.include_main.eq('t')]
mapping['clean']=mapping.icd_code.str.replace('.','',regex=False).str.strip().str.upper()
cardiac=set(mapping.loc[mapping.shock_class.eq('CARDIOGENIC_SHOCK'),'clean'])
bleed=set(mapping.loc[mapping.shock_class.eq('HEMORRHAGIC_OR_HYPOVOLEMIC_SHOCK'),'clean'])
later_admissions={}
for pid,group in cases.groupby('PatientID',sort=False):
    points=sorted(float(r.OffsetAfterFirstAdmission+r.ICUOffset) for r in group.itertuples(index=False) if finite(r.OffsetAfterFirstAdmission) and finite(r.ICUOffset))
    later_admissions[pid]=points

def case_info(case):
    E=float(case.TimeOfStay); D=float(case.OffsetOfDeath) if finite(case.OffsetOfDeath) else np.nan
    state=case.DischargeState; ht=case.HospitalDischargeType
    day=float(case.HospitalDischargeDay) if finite(case.HospitalDischargeDay) else np.nan
    hlo=day*86400 if finite(day) else np.nan; hhi=(day+1)*86400 if finite(day) else np.nan
    hosp_earlier=finite(hhi) and hhi<=E
    alive_conflict=finite(D) and D<=E and state==2202
    dead_conflict=finite(D) and D>E and state==2215
    hospital_alive_conflict=finite(D) and finite(hlo) and D<=hlo and ht in HOSP_ALIVE
    unknown_death=not finite(D) and (state==2215 or ht in HOSP_DEAD)
    terminal=finite(D) and D<=E
    hospital_normal=ht in HOSP_ALIVE and finite(hhi) and hhi>E
    normal=(state==2202 or hospital_normal) and not terminal and not unknown_death
    closure='DEATH' if terminal else ('NORMAL_CARE_END' if normal else 'CARE_END_UNCLASSIFIED')
    full=not any([hosp_earlier,alive_conflict,dead_conflict,hospital_alive_conflict]) and (terminal or normal)
    return {'E':E,'D':D,'hospital_lo':hlo,'hospital_hi':hhi,'hospital_type_alive':ht in HOSP_ALIVE,
            'hospital_earlier_conflict':hosp_earlier,'alive_death_conflict':alive_conflict,
            'dead_discharge_late_death_conflict':dead_conflict,'hospital_alive_death_conflict':hospital_alive_conflict,
            'known_death_time_missing':unknown_death,'closure':closure,'full_coverage':full}

def alive_evidence(case,ci,t):
    if finite(ci['D']):return ci['D']>t,'EXACT_LATER_DEATH'
    if ci['known_death_time_missing']:return False,'KNOWN_DEATH_TIME_UNRESOLVED'
    if case.DischargeState==2202 and ci['E']>=t:return True,'ALIVE_CASE_DISCHARGE'
    if ci['hospital_type_alive'] and finite(ci['hospital_lo']) and ci['hospital_lo']>=t:return True,'ALIVE_HOSPITAL_DAY_LOWER_BOUND'
    base=float(case.OffsetAfterFirstAdmission) if finite(case.OffsetAfterFirstAdmission) else np.nan
    if finite(base) and any(x>=base+t for x in later_admissions.get(case.PatientID,[])):
        return True,'DOCUMENTED_LATER_ADMISSION'
    return False,'NO_INDIVIDUAL_SURVIVAL_PROOF'

def outcome(row,case,ci):
    t=row['t0'];q=t+Q;h=t+H;E=ci['E'];D=ci['D'];r=row['next_raw_positive_start']
    validR=finite(r) and q<=r<=h and r<E and (not finite(D) or r<D)
    coverage_end=row['coverage_end']
    r_covered=validR and finite(coverage_end) and r<coverage_end
    # A known earlier restart is retained, even if later follow-up is incomplete.
    firstD=finite(D) and q<D<=h and D<=E and D<=coverage_end
    normal_end=ci['closure']=='NORMAL_CARE_END' and ci['full_coverage']
    observed_h=finite(coverage_end) and coverage_end>=h and E>=h
    binary=1 if r_covered else (0 if firstD or normal_end or observed_h else np.nan)
    if finite(D):death=int(t<D<=h);death_evidence='EXACT_DEATH_OFFSET'
    else:
        av,why=alive_evidence(case,ci,h)
        death=0 if av else np.nan;death_evidence=why
    comp=1 if binary==1 or death==1 else (0 if binary==0 and death==0 else np.nan)
    if r_covered:
        cif=True;st='RESTART';tm=r
    elif firstD:
        cif=True;st='DEATH';tm=D
    elif observed_h:
        cif=True;st='ADMIN_HORIZON';tm=h
    elif normal_end:
        cif=True;st='NORMAL_CARE_END';tm=E
    else:
        # A verified continuous record-level coverage end can be censored, but
        # neither an unclassified end nor absent rows establishes interruption.
        cif=finite(coverage_end) and coverage_end>=q
        st=('CARE_END_UNCLASSIFIED' if coverage_end==E else 'COVERAGE_END_UNCLASSIFIED') if cif else 'UNASSESSABLE_COVERAGE_END'
        tm=coverage_end if cif else np.nan
    ends=[(h,3,'ADMIN_HORIZON')]
    if finite(coverage_end):ends.append((coverage_end,2,'CARE_END_UNCLASSIFIED' if coverage_end==E else 'COVERAGE_END_UNCLASSIFIED'))
    if normal_end:ends.append((E,1,'NORMAL_CARE_END'))
    if firstD:ends.append((D,0,'DEATH'))
    end_time,_,end_reason=min(ends)
    row.update(restart_binary=binary,death_72h_binary=death,composite_restart_or_death=comp,
        death_outcome_evidence=death_evidence,cif_evaluable=cif,cif_state=st,
        cif_time_from_t0_hours=(tm-t)/3600 if finite(tm) else np.nan,
        valid_restart_time=r if r_covered else np.nan,followup_end_reason=end_reason,
        observation_end=E,ascertained_followup_end=end_time,
        zero_postqualification_observation=E==q,
        raw_restart_equal_q=finite(r) and r==q,raw_restart_equal_death=finite(r) and finite(D) and r==D,
        raw_restart_equal_care_end=finite(r) and r==E,raw_restart_equal_horizon=finite(r) and r==h)

rows=[]
for i,(cid,m) in enumerate(med.groupby('CaseID',sort=True)):
    case=valid.loc[cid];ci=case_info(case);origin=float(case.ICUOffset);E=ci['E'];weight=float(case.WeightOnAdmission)/1000
    edges=defaultdict(list);bad=[];raw_positive=[]
    for r in m[m.DrugID.isin(FIVE)].itertuples(index=False):
        if not finite(r.Offset) or not finite(r.OffsetDrugEnd) or not finite(r.AmountPerMinute):
            bad.append((r.Offset,r.OffsetDrugEnd));continue
        if r.AmountPerMinute<=0:continue
        if r.OffsetDrugEnd<=r.Offset:bad.append((r.Offset,r.OffsetDrugEnd));continue
        start=max(float(r.Offset),origin);stop=float(r.OffsetDrugEnd)
        if stop<=start:continue
        raw_positive.append((start,stop))
        if start>=E:continue
        rate=float(r.AmountPerMinute)*FACTOR[int(r.DrugID)]*(1 if r.DrugID==1550 else 1e6/weight)
        edges[start].append((1,int(r.id),int(r.DrugID),rate))
        edges[stop].append((-1,int(r.id),int(r.DrugID),rate))
    times=sorted(edges);active={d:{} for d in FIVE};segments=[]
    for j,a in enumerate(times[:-1]):
        for direction,rid,drug,rate in edges[a]:
            if direction==1:active[drug][rid]=rate
            else:active[drug].pop(rid,None)
        b=times[j+1];used=[d for d in FIVE if active[d]]
        if not used or b<=a:continue
        amb=any(len(active[d])>1 for d in used)
        dose=np.nan if amb else math.floor(sum(next(iter(active[d].values())) for d in used)*1e4+.5)/1e4
        if finite(dose) and dose<=0:
            QC['positive_rate_segments_rounded_to_zero']+=1
            continue
        segments.append({'start':a,'end':b,'dose':dose,'ambiguous':amb})
    raw_positive=merge_ranges(raw_positive)
    QC['malformed_five_drug_source_rows']+=len(bad)
    islands=[]
    for s in segments:
        if not islands or s['start']>=islands[-1][-1]['end']+Q:islands.append([s])
        else:islands[-1].append(s)
    lact=labgroups.get(int(cid));code=str(case.ICD10Main).replace('.','').strip().upper()
    for k,island in enumerate(islands):
        start=island[0]['start'];t=island[-1]['end'];q=t+Q
        dur=sum(s['end']-s['start'] for s in island)/3600
        peak=np.nan if any(s['ambiguous'] for s in island) else max(s['dose'] for s in island)
        lv=np.nan
        if lact is not None:
            v=lact[(lact.Offset>=start-21600)&(lact.Offset<=min(start+21600,t))]
            if len(v):lv=float(v.LaboratoryValue.max())
        ino=m[m.DrugID.isin([1559,1560])&(m.Offset<t)&(m.OffsetDrugEnd>start)&(m.AmountPerMinute>0)]
        prbc=m[(m.DrugID==2046)&(m.Offset>=start-86400)&(m.Offset<=min(start+86400,t))&(m.Amount>0)]
        flags={'positive_exposure_under_6h':dur<6,'peak_threshold_not_documented':not finite(peak) or peak<.1,
               'lactate_threshold_not_documented':not finite(lv) or lv<=2,'known_ECMO':cid in ecmo,
               'cardiogenic_main_ICD':code in cardiac,'inotrope_during_episode':len(ino)>0,
               'hemorrhagic_main_ICD':code in bleed,'PRBC_ge1400_ml':prbc.Amount.sum()>=1400}
        clinical=not any(flags.values())
        next_nee=islands[k+1][0]['start'] if k+1<len(islands) else np.nan
        r=next((a for a,b in raw_positive if a>=t),np.nan)
        raw_positive_in_qualification=any(a<q and b>t for a,b in raw_positive)
        row={'candidate_id':len(rows)+1,'CaseID':int(cid),'PatientID':int(case.PatientID),'episode_number':k+1,
             'episode_start':start,'t0':t,'q':q,'horizon':t+H,'care_end':E,'death':ci['D'],'interruption':np.nan,
             'cumulative_positive_hours':dur,'peak_nee':peak,'onset_lactate':lv,
             'clinical_gate_pass':clinical,'clinical_failure_reasons':'|'.join(z for z,v in flags.items() if v),
             'next_raw_positive_start':r,'next_rounded_nee_episode_start':next_nee,
             'raw_positive_in_qualification':raw_positive_in_qualification,
             'closure_class':ci['closure'],'coverage_type':'UNCONFIRMED',
             'coverage_end':np.nan,'hospital_day_before_case_conflict':ci['hospital_earlier_conflict'],
             'alive_discharge_death_conflict':ci['alive_death_conflict'],
             'dead_discharge_late_death_conflict':ci['dead_discharge_late_death_conflict'],
             'hospital_alive_death_conflict':ci['hospital_alive_death_conflict'],
             'death_time_missing_with_dead_status':ci['known_death_time_missing'],
             'source_interval_issues':len(bad),'raw_terminal_not_clipped':True}
        av,avwhy=alive_evidence(case,ci,q);row['alive_q_evidence']=avwhy
        fail=[];unk=[]
        if finite(ci['D']) and ci['D']<=q:
            if any(ci[z] for z in ['alive_death_conflict','dead_discharge_late_death_conflict','hospital_alive_death_conflict']):
                unk.append('DEATH_AT_OR_BEFORE_Q_WITH_SOURCE_CONFLICT')
            else:fail.append('DEATH_AT_OR_BEFORE_Q')
        if E<q:fail.append('CARE_END_BEFORE_Q')
        if raw_positive_in_qualification:fail.append('POSITIVE_DRUG_BEFORE_Q')
        if t>=E:fail.append('TERMINAL_SOURCE_NOT_BEFORE_CARE_END')
        if ci['full_coverage']:
            row['coverage_type']='FULL_MEDICATION_TABLE_AND_'+ci['closure'];row['coverage_end']=min(E,ci['D']) if finite(ci['D']) else E
        elif not any(ci[x] for x in ['alive_death_conflict','hospital_earlier_conflict','dead_discharge_late_death_conflict','hospital_alive_death_conflict']):
            ranges=[b for a,b in extra.get(int(cid),[]) if a<=t and b>=q]
            if ranges:
                row['coverage_type']='NON_TARGET_ADMINISTERED_INTERVAL_PLUS_INDEPENDENT_SURVIVAL'
                row['coverage_end']=min(max(ranges),E,ci['D'] if finite(ci['D']) else E)
        if not av:unk.append(avwhy)
        if ci['alive_death_conflict']:unk.append('ALIVE_DISCHARGE_DEATH_TIME_CONFLICT')
        if ci['dead_discharge_late_death_conflict']:unk.append('DEAD_DISCHARGE_DEATH_AFTER_CASE_CONFLICT')
        if ci['hospital_alive_death_conflict']:unk.append('ALIVE_HOSPITAL_DAY_AFTER_DEATH_CONFLICT')
        if ci['hospital_earlier_conflict']:unk.append('HOSPITAL_END_DAY_PRECEDES_CASE_CLOSE')
        if bad:unk.append('MALFORMED_SOURCE_INTERVALS')
        if not finite(row['coverage_end']) or row['coverage_end']<q:unk.append('QUALIFICATION_DRUG_COVERAGE_UNCONFIRMED')
        row['chronology_status']='FAIL' if fail else ('UNKNOWN' if unk else 'PASS')
        row['chronology_reasons']='|'.join(fail if fail else unk)
        row['eligible_candidate']=clinical and row['chronology_status']=='PASS'
        row['first_identity_blocking']=clinical and row['chronology_status']=='UNKNOWN'
        outcome(row,case,ci) if row['eligible_candidate'] else None
        rows.append(row)
    if i%2500==0:log(f'Cases processed={i}')

allc=pd.DataFrame(rows).sort_values(['CaseID','t0','episode_start']).reset_index(drop=True)
allc.to_csv(WORK/'all_source_candidate_episodes.csv.gz',index=False)
chosen=[];identity=[]
for cid,g in allc.groupby('CaseID',sort=True):
    possible=g[g.clinical_gate_pass&g.chronology_status.isin(['PASS','UNKNOWN'])]
    if not len(possible):identity.append({'CaseID':cid,'selection_status':'NO_ELIGIBLE_CANDIDATE'});continue
    z=possible.iloc[0]
    if z.chronology_status=='UNKNOWN':identity.append({'CaseID':cid,'selection_status':'FIRST_IDENTITY_UNKNOWN','first_candidate_id':z.candidate_id})
    else:
        chosen.append(z.to_dict());identity.append({'CaseID':cid,'selection_status':'FIRST_ELIGIBLE_SELECTED','first_candidate_id':z.candidate_id})
sel=pd.DataFrame(chosen);ident=pd.DataFrame(identity)
sel.to_csv(WORK/'selected_cohort.csv.gz',index=False)
ident.to_csv(WORK/'first_identity_by_case.csv.gz',index=False)
