"""Convert audited minimal-cohort states to CIF input; no outcome adjudication."""
from pathlib import Path
from datetime import datetime, timedelta
from collections import Counter
import csv
import hashlib
import json

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
R = private_output('mimic_clinical')
P = R.parents[1]
A = R / 'aggregate'
W = R / 'restricted_work'
A.mkdir(exist_ok=True)
W.mkdir(exist_ok=True)
assert not (A / 'input_qc.json').exists(), 'Completed preparation is protected'

paths={'audited_selected_cohort':release_path('mimic_selected_cohort'), 'raw_observed_predictors':release_path('mimic_observed_input')}

def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def write(p, rows):
    with p.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def group(x):
    return '<0.2' if x < .2 else ('0.2-<0.5' if x < .5 else '>=0.5')

before = {k: sha(p) for k, p in paths.items()}
selected = read(paths['audited_selected_cohort'])
raw = read(paths['raw_observed_predictors'])
assert len(selected) == 1958 and len(raw) == 1958
raw_by_stay = {r['stay_id']: r for r in raw}
assert len(raw_by_stay) == 1958
assert {r['stay_id'] for r in selected} == set(raw_by_stay)
allowed = {'RESTART', 'DEATH', 'NORMAL_CARE_END', 'ADMIN_HORIZON'}
out = []
for r in selected:
    s = raw_by_stay[r['stay_id']]
    t0 = datetime.fromisoformat(r['t0'])
    assert datetime.fromisoformat(s['t0']) == t0
    assert s['subject_id'] == r['subject_id'] and s['hadm_id'] == r['hadm_id']
    assert r['time_status'] == 'PASS' and r['cif_status'] in allowed
    assert r['qualification_time'] and r['horizon'] and r['cif_time']
    assert datetime.fromisoformat(r['qualification_time']) == t0 + timedelta(hours=4)
    assert datetime.fromisoformat(r['horizon']) == t0 + timedelta(hours=72)
    tm = (datetime.fromisoformat(r['cif_time']) - t0).total_seconds() / 3600
    assert 4 <= tm <= 72
    peak = float(s['predictor_episode_peak_nee'])
    assert abs(peak - float(r['peak_nee'])) < 1e-12
    st = r['cif_status']
    main = {'RESTART': 1, 'DEATH': 2}.get(st, 0)
    s20 = {'RESTART': 1, 'DEATH': 2, 'NORMAL_CARE_END': 3}.get(st, 0)
    if st == 'RESTART':
        assert r['restart_72h'] == '1'
        assert datetime.fromisoformat(r['cif_time']) < datetime.fromisoformat(r['care_end'])
        assert not r['deathtime'] or datetime.fromisoformat(r['cif_time']) < datetime.fromisoformat(r['deathtime'])
    if s20 == 3:
        assert r['care_end_type'] == 'NORMAL_CARE_END'
    out.append(dict(cohort='MIMIC', version='minimal_1958', subject_id=r['subject_id'],
                    stay_id=r['stay_id'], time_hours=tm, status_main=main, status_s20=s20,
                    state=st, peak_nee=peak, dose_group=group(peak), entry_hours=4))
assert len({r['subject_id'] for r in out}) == 1856
assert Counter(r['state'] for r in out) == {'RESTART': 680, 'DEATH': 30, 'NORMAL_CARE_END': 681, 'ADMIN_HORIZON': 567}
write(W / 'clinical_course_inputs.csv', out)
counts = []
for g in ['Overall', '<0.2', '0.2-<0.5', '>=0.5']:
    z = out if g == 'Overall' else [r for r in out if r['dose_group'] == g]
    counts.append(dict(cohort='MIMIC', version='minimal_1958', dose_group=g, n=len(z),
                       patients=len({r['subject_id'] for r in z}),
                       first_restart=sum(r['status_main'] == 1 for r in z),
                       first_death=sum(r['status_main'] == 2 for r in z),
                       normal_end_competitor=sum(r['status_s20'] == 3 for r in z),
                       primary_censored=sum(r['status_main'] == 0 for r in z),
                       s20_censored=sum(r['status_s20'] == 0 for r in z),
                       immediate_restart_at_q=sum(r['status_main'] == 1 and r['time_hours'] == 4 for r in z)))
write(A / 'clinical_course_event_counts.csv', counts)
assert all(sha(paths[k]) == h for k, h in before.items())
write(A / 'source_manifest.csv', [dict(role=k, path=str(paths[k]), sha256=h, unchanged=True) for k, h in before.items()])
qc = dict(PASS=True, events=len(out), patients=1856, raw_peak_max_difference=0,
          source_statuses_reused_not_readjudicated=True, unknown_promoted_to_normal=False,
          input_only_no_model_results_read=True, external_cohorts_read=False,
          entry_hours=4, horizon_hours=72, immediate_restart_at_q=counts[0]['immediate_restart_at_q'],
          before_q_rows=0, source_hashes_unchanged=True,
          input_sha256=sha(W / 'clinical_course_inputs.csv'))
(A / 'input_qc.json').write_text(json.dumps(qc, indent=2), encoding='utf-8')
print(json.dumps(qc, indent=2))
