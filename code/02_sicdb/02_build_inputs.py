"""Extract all predictors for a fixed final cohort using the source-derived SQLite store.

No historical cohort, feature cache, fitted object or row prediction is read.
The source window functions are retained in feature_functions.py.
"""
from pathlib import Path
import importlib.util
import sqlite3
import sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output

PRED=['predictor_age','predictor_female','predictor_episode_peak_nee','predictor_last_positive_nee',
      'predictor_positive_nee_hours','predictor_multivaso_peak','predictor_sofa_24h',
      'predictor_mechanical_ventilation_t0','predictor_rrt_t0','predictor_mean_map_4h',
      'predictor_map_source_invasive','predictor_mean_heart_rate_4h','predictor_urine_output_24h_ml']

def main():
    q=pd.read_csv(release_path('sicdb_selected_cohort')).set_index('CaseID',drop=False)
    if q.index.has_duplicates:raise ValueError('One fixed event per care case is required.')
    output=private_output('sicdb_inputs')
    if (output/'clinical_inputs_all.csv').exists():raise FileExistsError('Private output already exists.')
    spec=importlib.util.spec_from_file_location('features',Path(__file__).with_name('feature_functions.py'))
    features=importlib.util.module_from_spec(spec);spec.loader.exec_module(features)
    connection=sqlite3.connect(release_path('sicdb_sqlite').resolve().as_uri()+'?mode=ro',uri=True)
    connection.execute('PRAGMA query_only=ON')
    cases=pd.read_sql_query('SELECT CaseID,AgeOnAdmission,Sex,ICUOffset,TimeOfStay,WeightOnAdmission,HoursOfCRRT FROM cases',connection).set_index('CaseID')
    features.c=connection;features.cases=cases
    try:
        extracted=[]
        for cid,z in q.iterrows():
            values=features.episode_features(cid,z.episode_start,z.t0)|features.window_features(cid,z.t0)
            values['CaseID']=cid;extracted.append(values)
        fresh=pd.DataFrame(extracted).set_index('CaseID')
        for k in fresh.columns:q[k]=fresh[k]
    finally:
        connection.close()
    q['case_icu_start']=cases.loc[q.index,'ICUOffset']
    q['HoursOfCRRT']=cases.loc[q.index,'HoursOfCRRT']
    q['admission_weight_kg_binned']=cases.loc[q.index,'WeightOnAdmission']/1000
    q['case_end']=q.care_end
    q['model_row_id']=np.arange(1,len(q)+1)
    q['subject_id']=q.PatientID;q['stay_id']=q.CaseID;q['hadm_id']=q.CaseID
    q['d2v2_event_id']=q.candidate_id
    q['outcome_restart_72h']=q.restart_binary;q['evaluable']=q.restart_binary.notna()
    q['death_evaluable']=q.death_72h_binary.notna()
    q['composite_evaluable']=q.composite_restart_or_death.notna()
    q['predictor_rrt_t0_recorded']=q.recorded_crrt_active_pre6h
    q['rrt_clinical_state_unknown']=q.rrt_preonly_proxy.isna()
    q['first_status_main']=q.cif_state.map({'RESTART':1,'DEATH':2}).fillna(0).astype(int)
    q['first_time_hours']=q.cif_time_from_t0_hours
    q['first_time_since_q_hours']=q.cif_time_from_t0_hours-4
    q['first_status_case_competing']=np.where(q.cif_state.eq('NORMAL_CARE_END'),3,q.first_status_main)
    q['hours_from_first_ICU_IMC_entry']=(q.t0-q.case_icu_start)/3600
    q['death_hours']=(q.death-q.t0)/3600;q['case_end_hours']=(q.care_end-q.t0)/3600
    q['restart_hours']=(q.valid_restart_time-q.t0)/3600
    q.to_csv(output/'clinical_inputs_all.csv',index=False)
    required=[PRED[0],PRED[1],PRED[2],PRED[4],PRED[7],'predictor_rrt_t0_recorded']
    subsets={'association_inputs':required,'sofa_inputs':required+[PRED[6]],'scorable_inputs':PRED[:9]}
    for name,columns in subsets.items():
        data=q[q.evaluable&q[columns].notna().all(axis=1)].copy()
        if name!='scorable_inputs':data['predictor_rrt_t0']=data.predictor_rrt_t0_recorded
        if name=='association_inputs':data=data.drop(columns=['predictor_sofa_24h'])
        data.to_csv(output/(name+'.csv'),index=False)
    print('Private predictor inputs written. Outcomes were read, not reassigned; no model was fitted.')

if __name__=='__main__':
    main()
