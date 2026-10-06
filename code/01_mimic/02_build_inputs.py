"""Extract predictors for the fixed final cohort, without historical feature caches.

This does not select a cohort, assign outcomes, fit models or modify source tables.
PostgreSQL connection fields come from MIMIC_PG_* environment variables.
"""
from pathlib import Path
from datetime import datetime
import csv
import os
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output

def main():
    import psycopg
    with release_path('mimic_selected_cohort').open(encoding='utf-8-sig',newline='') as f:
        cohort = sorted(csv.DictReader(f),key=lambda r:(int(r['stay_id']),r['t0']))
    if not cohort or len({r['stay_id'] for r in cohort}) != len(cohort):
        raise ValueError('Expected one fixed cohort event per care stay.')
    requested=[]
    for i,r in enumerate(cohort,1):
        requested.append((i,-i,int(r['subject_id']),int(r['hadm_id']),int(r['stay_id']),
                          datetime.fromisoformat(r['episode_start']),datetime.fromisoformat(r['t0']),
                          float(r['positive_hours']),float(r['peak_nee']),int(r['restart_72h']),
                          int(r['death_72h']) if r['death_72h'] else None))
    columns='model_row_id,d2v2_event_id,subject_id,hadm_id,stay_id,episode_start,t0,positive_hours,peak_nee,outcome_restart_72h,secondary_death_72h'
    value_pattern='(%s::integer,%s::integer,%s::integer,%s::integer,%s::integer,%s::timestamp,%s::timestamp,%s::double precision,%s::double precision,%s::integer,%s::integer)'
    query='WITH requested('+columns+') AS (VALUES '+','.join([value_pattern]*len(requested))+'),\n'
    query+=Path('code/01_mimic/feature_query.sql').read_text(encoding='utf-8')
    config={k:os.environ['MIMIC_PG_'+k.upper()] for k in ('host','port','dbname','user','password')
            if os.environ.get('MIMIC_PG_'+k.upper())}
    output=private_output('mimic_inputs')
    raw_path=output/'raw_observed_modeling_table.csv'
    main_path=output/'main_modeling_table.csv'
    if raw_path.exists() or main_path.exists():
        raise FileExistsError('Private input outputs already exist; choose an empty destination.')
    with psycopg.connect(**config) as connection:
        connection.read_only=True
        with connection.cursor() as cursor:
            cursor.execute(query,[v for row in requested for v in row])
            names=[d.name for d in cursor.description]
            rows=[dict(zip(names,r)) for r in cursor.fetchall()]
        connection.rollback()
    if len(rows)!=len(cohort):
        raise ValueError('Feature join changed the number of cohort events.')
    # Ledger fields are not predictor inputs. Their original extraction logic stays in the SQL.
    names=[k for k in names if not k.startswith('ledger_')]
    with release_path('frozen_winsor').open(encoding='utf-8-sig',newline='') as f:
        clipping=list(csv.DictReader(f))
    clipped=[]
    for r in rows:
        z={k:r[k] for k in names}
        group=r['validation_anchor_year_group']
        if group in ('2008 - 2010','2011 - 2013','2014 - 2016'):
            role='DEVELOPMENT_2008_2016'
        elif group in ('2017 - 2019','2020 - 2022'):
            role='TEMPORAL_VALIDATION_2017_2022'
        else:
            raise ValueError('Unrecognised temporal group.')
        z['temporal_validation_role']=role
        for c in clipping:
            k=c['predictor']
            if z.get(k) is not None:
                z[k]=min(max(float(z[k]),float(c['lower_p005'])),float(c['upper_p995']))
        clipped.append(z)
    main_names=[k for k in names if k!='predictor_peak_lactate_12h']
    for path,data,fields in [(raw_path,rows,names),(main_path,clipped,main_names+['temporal_validation_role'])]:
        with path.open('x',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
            writer.writeheader();writer.writerows(data)
    print('Private model inputs written. No model was fitted.')

if __name__=='__main__':
    main()
