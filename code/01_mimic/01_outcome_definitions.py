"""Final-analysis functions; arguments and globals are documented by the source signatures."""
from pathlib import Path
from datetime import datetime, timedelta, time
from collections import defaultdict, Counter
import csv, json, hashlib, argparse

Q=timedelta(hours=4); H=timedelta(hours=72); DAY=timedelta(days=1)
ALIVE_DESTINATIONS={'HOME','HOME HEALTH CARE','HOSPICE','SKILLED NURSING FACILITY','REHAB','ACUTE HOSPITAL','CHRONIC/LONG TERM ACUTE CARE','ASSISTED LIVING','HEALTHCARE FACILITY','AGAINST ADVICE'}

def dt(x):return datetime.fromisoformat(x) if x else None

def boolean(x):return x in ('t','True','true','1')

def normalise(r):
    r=dict(r)
    for k in ('episode_start','t0','source_last_end','next_positive','intime','outtime','admittime','dischtime','deathtime','dod'):r[k]=dt(r[k])
    for k in ('subject_id','hadm_id','stay_id','episode','hospital_expire_flag','n_sepsis'):r[k]=int(r[k])
    for k in ('positive_hours','peak_nee','prbc_ml'):r[k]=float(r[k])
    r['peak_lactate']=float(r['peak_lactate']) if r['peak_lactate'] else None
    r['inotrope']=boolean(r['inotrope'])
    return r

def terminal_conflicts(r):
    d=r['deathtime'];day=r['dod'];loc=r['discharge_location'];expire=r['hospital_expire_flag']
    reasons=[]
    if expire==0 and loc=='DIED':reasons.append('ALIVE_FLAG_DEAD_DESTINATION')
    if expire==1 and loc in ALIVE_DESTINATIONS:reasons.append('DEAD_FLAG_ALIVE_DESTINATION')
    if d and d>r['dischtime']:reasons.append('DEATH_AFTER_HOSPITAL_CLOSE')
    if d and d<r['admittime']:reasons.append('DEATH_BEFORE_ADMISSION')
    if d and d<=r['dischtime'] and expire==0:reasons.append('EXACT_INHOSPITAL_DEATH_ALIVE_FLAG')
    if not d and day and day>r['dischtime'] and (expire==1 or loc=='DIED'):
        reasons.append('DEAD_AT_CLOSE_BUT_LATER_DEATH_DATE')
    return reasons

def qualification(r):
    t=r['t0']; q=t+Q; d=r['deathtime']; day=r['dod']; e=min(r['outtime'],r['dischtime'])
    if e<t:return 'FAIL','T0_AFTER_CARE_END'
    if e<q:return 'FAIL','CARE_END_BEFORE_Q'
    if r['source_last_end']!=t:return 'FAIL','ADMINISTRATIVELY_TRUNCATED_END'
    if r['next_positive'] and r['next_positive']<q:return 'FAIL','RESTART_BEFORE_Q'
    if r.get('rounded_zero_drug_in_qualification',False):return 'FAIL','POSITIVE_DRUG_ROUNDED_ZERO_WITHIN_Q'
    if r['intime']>t or r['admittime']>t:return 'UNKNOWN','INVALID_ENTRY_ORDER'
    if terminal_conflicts(r):return 'UNKNOWN','UNRESOLVED_TERMINAL_SOURCE_CONFLICT'
    if d and day and d.date()!=day.date():
        if max(d,day+DAY)<=q:return 'FAIL','BOTH_DEATH_SOURCES_BEFORE_Q'
        if min(d,day)<=q:return 'UNKNOWN','DEATH_DATE_SOURCE_CONFLICT'
    if d and d<=t:return 'FAIL','DEATH_AT_OR_BEFORE_T0'
    if d and d<=q:return 'FAIL','DEATH_BEFORE_OR_AT_Q'
    if not d and day:
        if day+DAY<=t:return 'FAIL','DATE_DEATH_BEFORE_T0'
        if day+DAY<=q:return 'FAIL','DATE_DEATH_BEFORE_Q'
        if day<=q:return 'UNKNOWN','DATE_DEATH_OVERLAPS_QUALIFICATION'
    if not d and not day and (r['hospital_expire_flag']==1 or r['discharge_location']=='DIED'):
        return 'UNKNOWN','KNOWN_DEATH_NO_TIME'
    # Whole-source extraction, valid stay/admission and terminal medication evidence.
    # No timestamps are imputed; no arbitrary measurement-density rule is introduced.
    return 'PASS','QUALIFIED_RECORD_COVERAGE'

def evaluate(r,lastrecord):
    t=r['t0'];q=t+Q;h=t+H;d=r['deathtime'];day=r['dod'];e=min(r['outtime'],r['dischtime']);nr=r.get('effective_next_positive',r['next_positive'])
    source_conflict=bool(d and day and d.date()!=day.date())
    relevant_source_conflict=bool(source_conflict and min(d,day)<=min(e,h))
    # Classify known normal care ending independently from mortality follow-up.
    if d and d<=e:ending='DEATH'
    elif d and d>e:ending='NORMAL_CARE_END'
    elif day and day>e:ending='NORMAL_CARE_END'
    elif day and day<=e<day+DAY:ending='CARE_END_UNCLASSIFIED'
    elif r['hospital_expire_flag']==0 and r['discharge_location']!='DIED':ending='NORMAL_CARE_END'
    elif r['hospital_expire_flag']==1 and not d:ending='CARE_END_UNCLASSIFIED'
    else:ending='NORMAL_CARE_END' if e==r['outtime'] and e<r['dischtime'] else 'CARE_END_UNCLASSIFIED'
    stop=min([h,e]+([d] if d else []))
    coverage='FULL_SOURCE_STRUCTURED_CARE_NO_KNOWN_INTERRUPTION'
    # Whether no observed restart is a genuine negative within the defined care scope.
    date_at_stop=bool(not d and day and day<=min(e,h))
    plausible_r=bool(nr and q<=nr<=h and nr<e)
    restart=None; restart_reason='UNCLASSIFIED_END_WITHOUT_RESTART'
    valid_r=bool(plausible_r and (not d or nr<d) and (not day or (d and not source_conflict) or nr<day))
    ambiguous_r=bool(plausible_r and not d and day and day<=nr<day+DAY)
    if valid_r:restart=1;restart_reason='VALID_RESTART'
    elif ambiguous_r or relevant_source_conflict:restart=None;restart_reason='RESTART_DEATH_ORDER_UNKNOWN'
    elif (d and q<d<=min(e,h)) or ending=='NORMAL_CARE_END' or (e>=h and not date_at_stop):
        restart=0;restart_reason='NO_RESTART_WITH_DEFINITE_ENDPOINT'
    elif date_at_stop and (not plausible_r or nr>=day+DAY):
        # Entire relevant medication table is observed; no positive interval before care end.
        # Death ordering may remain unavailable to CIF even when binary restart is known.
        restart=0;restart_reason='NO_POSITIVE_RECORD_BEFORE_DATE_DEATH_OR_CARE_END'
    # Independent 72h all-cause mortality, never infer survival solely from missing death.
    death=None
    if d:death=int(t<d<=h)
    elif day:
        if day>t and day+DAY<=h:death=1
        elif day>h:death=0
    elif lastrecord and lastrecord>=h:death=0
    elif e>=h and r['hospital_expire_flag']==0:death=0
    if source_conflict:
        date_status=1 if day>t and day+DAY<=h else (0 if day>h else None)
        if date_status is None or death!=date_status:death=None
    composite=1 if restart==1 or death==1 else (0 if restart==0 and death==0 else None)
    # Risk-set entry q; no estimated CIF in this task.
    cif='UNKNOWN';cif_time=None
    if valid_r:cif='RESTART';cif_time=nr
    elif d and q<d<=h and d<=e and not relevant_source_conflict:cif='DEATH';cif_time=d
    elif relevant_source_conflict:cif='UNKNOWN_DEATH_SOURCE_ORDER'
    elif not d and day and day<=min(e,h):cif='UNKNOWN_DATE_OR_ORDER'
    elif e>=h:cif='ADMIN_HORIZON';cif_time=h
    elif ending=='NORMAL_CARE_END':cif='NORMAL_CARE_END';cif_time=e
    else:cif='CARE_END_UNCLASSIFIED'
    if ambiguous_r:cif='UNKNOWN_DATE_OR_ORDER';cif_time=None
    return dict(qualification_time=q,horizon=h,observation_end=e,ascertained_followup_end=None if relevant_source_conflict or date_at_stop else stop,
                followup_end_reason='UNKNOWN_DEATH_TIME' if relevant_source_conflict or date_at_stop else ('DEATH' if d and stop==d else ('ADMIN_HORIZON' if stop==h else ending)),
                care_end=e,care_end_type=ending,interruption_time=None,interruption_status='NO_EXPLICIT_FEED_UPTIME_FIELD',coverage_basis=coverage,
                restart_72h=restart,restart_reason=restart_reason,death_72h=death,composite_72h=composite,cif_status=cif,cif_time=cif_time,
                zero_postqualification_observation=int(e==q),recorded_restart_at_death=int(bool(nr and d and nr==d)),recorded_restart_at_care_end=int(bool(nr and nr==e)),
                rejected_postdeath_restart=int(bool(plausible_r and d and nr>=d)),rejected_at_or_after_care_end_restart=int(bool(nr and q<=nr<=h and nr>=e)))
