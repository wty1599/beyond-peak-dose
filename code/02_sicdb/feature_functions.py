"""Original episode and pre-withdrawal feature definitions; caller supplies c and cases."""
from collections import defaultdict
import math
import numpy as np
import pandas as pd
F={1562:1.,1502:1.,1618:.01,1593:.1,1550:2.5}
VENT=[712,713,715,717,718,2019,2278,2284,3040]
def episode_features(cid,start,t):
 ca=cases.loc[cid];weight=float(ca.WeightOnAdmission)/1000
 meds=pd.read_sql_query('SELECT id,DrugID,Offset,OffsetDrugEnd,AmountPerMinute FROM medication WHERE CaseID=? AND DrugID IN (1562,1502,1618,1593,1550) AND Offset<? AND OffsetDrugEnd>?',c,params=(int(cid),float(t),float(start)))
 edges=defaultdict(list)
 for m in meds.itertuples(index=False):
  if not(m.AmountPerMinute>0 and m.OffsetDrugEnd>m.Offset):continue
  left=max(float(m.Offset),float(ca.ICUOffset),start);right=min(float(m.OffsetDrugEnd),t)
  if right<=left:continue
  rate=float(m.AmountPerMinute)*F[m.DrugID]
  if m.DrugID!=1550:rate=rate*1e6/weight
  edges[left].append((1,m.id,m.DrugID,rate));edges[right].append((-1,m.id,m.DrugID,rate))
 ts=sorted(edges);active={d:{} for d in F};segments=[]
 for i,left in enumerate(ts[:-1]):
  for direction,rid,drug,rate in edges[left]:
   if direction==1:active[drug][rid]=rate
   else:active[drug].pop(rid,None)
  used=[d for d in F if active[d]]
  if not used:continue
  assert all(len(active[d])==1 for d in used),'Unresolved same-drug overlap'
  dose=math.floor(sum(next(iter(active[d].values())) for d in used)*10000+.5)/10000
  if dose>0:segments.append((left,ts[i+1],dose,len(used)))
 x=np.asarray(segments)
 assert len(x) and x[0,0]==start and x[-1,1]==t
 peak=x[:,2].max()
 return {'predictor_age':float(ca.AgeOnAdmission),'predictor_female':{735:0,736:1}.get(ca.Sex,np.nan),
  'predictor_episode_peak_nee':peak,'predictor_last_positive_nee':x[-1,2],
  'predictor_positive_nee_hours':np.sum(x[:,1]-x[:,0])/3600,
  'predictor_multivaso_peak':float(x[x[:,2]==peak,3].max()>=2)}

def window_features(cid,t):
 sig=pd.read_sql_query('SELECT DataID,Offset,Val,cnt FROM signals WHERE CaseID=? AND Offset>=? AND Offset<=?',c,params=(int(cid),float(t)-172800,float(t)))
 sig=sig[np.isfinite(sig.Val)&(sig.cnt>0)].copy();sig['end']=sig.Offset+3600
 pre4=sig[(sig.Offset>=t-14400)&(sig.end<=t)]
 pre24=sig[(sig.Offset>=t-86400)&(sig.end<=t)]
 inv=pre4[(pre4.DataID==703)&(pre4.Val>0)&(pre4.Val<300)]
 ni=pre4[(pre4.DataID==706)&(pre4.Val>0)&(pre4.Val<300)]
 maps=inv if len(inv) else ni
 hr=pre4[pre4.DataID.isin([707,708,724])&(pre4.Val>0)&(pre4.Val<300)]
 hr=hr.assign(priority=hr.DataID.map({707:0,708:1,724:2})).sort_values(['Offset','priority']).drop_duplicates('Offset')
 urine=pre24[(pre24.DataID==725)&(pre24.Val>=0)]
 sofa=sig[(sig.DataID==3139)&(sig.Offset+86400<=t)&(sig.Offset+86400>t-86400)&(sig.Val>=0)&(sig.Val<=24)].sort_values('Offset')
 airway=c.execute('SELECT EXISTS(SELECT 1 FROM data_range WHERE CaseID=? AND DataID IN (720,3041) AND Offset<=? AND OffsetEnd>?)',(int(cid),float(t),float(t))).fetchone()[0]
 vent=pre4[pre4.DataID.isin(VENT)&(pre4.Val>0)]
 observed=pre4[pre4.DataID.isin([703,706,707,708,724])&(pre4.Val>0)]
 mv=1 if airway and len(vent) else (0 if not airway and not len(vent) and len(observed) else np.nan)
 rr=sig[sig.DataID.isin([723,730,731,732,2022])&(sig.Offset>t-21600)&(sig.end<=t)&(sig.Val>=0)]
 blood=rr[rr.DataID==723].sort_values('Offset');state=np.nan;end=np.nan
 if len(blood):
  last=blood.iloc[-1];later=rr[rr.Offset>=last.Offset];state=float(last.Val>0 or (later.Val>0).any());end=float(later.end.max())
 elif (rr.Val>0).any():state=1.;end=float(rr.loc[rr.Val>0,'end'].max())
 hist=c.execute('SELECT EXISTS(SELECT 1 FROM signals WHERE CaseID=? AND DataID IN (723,730,731,732) AND Val>0)',(int(cid),)).fetchone()[0]
 legacy=state
 if pd.isna(legacy) and cases.loc[cid,'HoursOfCRRT']==0 and not hist:legacy=0.
 return {'predictor_sofa_24h':float(sofa.iloc[-1].Val) if len(sofa) else np.nan,
  'sofa_label_offset':float(sofa.iloc[-1].Offset) if len(sofa) else np.nan,
  'sofa_assumed_window_end':float(sofa.iloc[-1].Offset+86400) if len(sofa) else np.nan,
  'predictor_mechanical_ventilation_t0':mv,'predictor_rrt_t0':legacy,
  'predictor_mean_map_4h':float(maps.Val.mean()) if len(maps) else np.nan,
  'predictor_map_source_invasive':1. if len(inv) else (0. if len(ni) else np.nan),
  'predictor_mean_heart_rate_4h':float(hr.Val.mean()) if len(hr) else np.nan,
  'predictor_urine_output_24h_ml':float(urine.Val.sum()) if len(urine) else np.nan,
  'rrt_preonly_proxy':state,'recorded_crrt_active_pre6h':0. if pd.isna(state) else state,
  'recorded_crrt_evidence':('no_activity_recorded_NOT_confirmed_no_RRT' if pd.isna(state) else ('activity_recorded' if state else 'recent_zero_anchor')),
  'rrt_max_used_bin_end':end}
