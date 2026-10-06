"""Final-analysis functions; arguments and globals are documented by the source signatures."""
from pathlib import Path
from collections import defaultdict, Counter
import argparse, hashlib, json, math, sqlite3, sys, time
import numpy as np
import pandas as pd

Q=14400; H=72*3600
HOSP_ALIVE={2026,3129,3131,3132,3133,3134}; HOSP_DEAD={2028,3130}
QC=Counter()

def finite(x): return x is not None and pd.notna(x) and np.isfinite(x)

def merge_ranges(ranges):
    out=[]
    for a,b in sorted(ranges):
        if not out or a>out[-1][1]:out.append([a,b])
        else:out[-1][1]=max(out[-1][1],b)
    return out

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
