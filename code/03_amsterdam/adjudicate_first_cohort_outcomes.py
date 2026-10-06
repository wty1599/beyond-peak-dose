"""Final database-specific algorithm. All clinical outputs are private and must not be committed."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output
import build_technical_episodes as technical
ROOT = Path(__file__).resolve().parent
STAGE = private_output('amsterdam_cohort')
ADMISSIONS = release_path('amsterdam_raw') / 'admissions.csv'
VASO = STAGE / 'vaso_source.csv.gz'
PROCESS = STAGE / 'process_source.csv.gz'
PROXY = STAGE / 'clinical_proxy_candidates_v1.csv.gz'
PROXY_AUDIT = STAGE / 'clinical_proxy_candidates_v1_audit.json'
OUT = STAGE / 'first_cohort_outcomes_v1.csv.gz'
PART = STAGE / 'first_cohort_outcomes_v1.csv.gz.incomplete'
AUDIT = STAGE / 'first_cohort_outcomes_v1_audit.json'
DAY = 86400000
HOUR = 3600000
SUPPORT_IDS = {'ventilation_process': '9328', 'cvvh_process': '12465', 'hemodialysis_process': '16363'}

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def integer(value: str | None) -> int | None:
    try:
        return int(value) if value not in (None, '') else None
    except (ValueError, TypeError):
        return None

def destination_type(value: str) -> str:
    if value == 'Overleden':
        return 'DEATH_TYPE'
    if value and integer(value) is not None:
        return 'RECORDED_NONDEATH_END'
    return 'UNKNOWN_END_TYPE'

def death_envelope(date_value: int | None) -> tuple[int, int] | None:
    return None if date_value is None else (date_value - DAY, date_value + DAY)

def check_inputs() -> None:
    for path in (ADMISSIONS, VASO, PROCESS, PROXY):
        if path.is_symlink() or not path.is_file():
            raise RuntimeError('Required private input missing')
    if not PROXY_AUDIT.is_file() or PROXY_AUDIT.is_symlink():
        raise RuntimeError('proxy_audit_missing')
    audit = json.loads(PROXY_AUDIT.read_text(encoding='utf-8'))
    if audit.get('no_outcomes_calculated') is not True or audit.get('technical_episode_rows') != 15533:
        raise RuntimeError('proxy_audit_mismatch')
    for target in (OUT, PART, AUDIT):
        if target.exists():
            raise RuntimeError(f'refuse_overwrite:{target.name}')

def read_admissions() -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]]]:
    by_id = {}
    by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    with ADMISSIONS.open('r', encoding='cp1252', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        needed = {'patientid', 'admissionid', 'admittedat', 'dischargedat', 'destination', 'dateofdeath'}
        if not needed.issubset(reader.fieldnames or []):
            raise RuntimeError('admissions_outcome_schema_mismatch')
        for row in reader:
            if row['admissionid'] in by_id:
                raise RuntimeError('duplicate_admission_id')
            record = {field: row[field] for field in needed}
            by_id[row['admissionid']] = record
            by_patient[row['patientid']].append(record)
    if len(by_id) != 23106:
        raise RuntimeError('admissions_row_count_mismatch')
    return (by_id, by_patient)

def read_gzip_by_admission(path: Path, expected_count: int) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        if 'admissionid' not in (reader.fieldnames or []):
            raise RuntimeError(f'missing_admissionid:{path.name}')
        for count, row in enumerate(reader, 1):
            if count > expected_count:
                raise RuntimeError(f'input_row_bound_exceeded:{path.name}')
            out[row['admissionid']].append(row)
    if count != expected_count:
        raise RuntimeError(f'input_row_count_mismatch:{path.name}')
    return out

def patient_death(rows: list[dict[str, str]]) -> tuple[tuple[int, int] | None, bool]:
    dates: set[int] = set()
    conflict = False
    for row in rows:
        raw = row['dateofdeath']
        if not raw:
            continue
        parsed = integer(raw)
        if parsed is None:
            conflict = True
        else:
            dates.add(parsed)
    if conflict or len(dates) > 1:
        return (None, True)
    return (death_envelope(next(iter(dates))) if dates else None, False)

def later_living_admission(rows: list[dict[str, str]], current_admission: str, instant: int) -> bool:
    return any((row['admissionid'] != current_admission and (later := integer(row['admittedat'])) is not None and (later > instant) for row in rows))

def source_conflict(admission: dict[str, str], patient_rows: list[dict[str, str]], death: tuple[int, int] | None, death_conflict: bool) -> bool:
    if death_conflict:
        return True
    c = integer(admission['dischargedat'])
    if c is None:
        return True
    end_type = destination_type(admission['destination'])
    if death is None:
        return end_type == 'DEATH_TYPE' and later_living_admission(patient_rows, admission['admissionid'], c)
    low, high = death
    admitted = integer(admission['admittedat'])
    if admitted is not None and high <= admitted:
        return True
    if end_type == 'RECORDED_NONDEATH_END' and high <= c:
        return True
    if end_type == 'DEATH_TYPE' and c < low:
        return True
    return any((row['admissionid'] != admission['admissionid'] and (later := integer(row['admittedat'])) is not None and (high <= later) for row in patient_rows))

def survival_record_after(admission: dict[str, str], patient_rows: list[dict[str, str]], instant: int) -> bool:
    c = integer(admission['dischargedat'])
    nondeath_end = destination_type(admission['destination']) == 'RECORDED_NONDEATH_END'
    return nondeath_end and c is not None and (c > instant) or later_living_admission(patient_rows, admission['admissionid'], instant)

def q_qualification(admission: dict[str, str], patient_rows: list[dict[str, str]], t0: int, death: tuple[int, int] | None, conflict: bool) -> str:
    q = t0 + 4 * HOUR
    c = integer(admission['dischargedat'])
    if c is None:
        return 'UNKNOWN_Q_CARE_BOUNDARY'
    if c < q:
        return 'FAIL_CARE_END_BEFORE_Q'
    if conflict:
        return 'UNKNOWN_DEATH_SOURCE_CONFLICT'
    if death is not None:
        low, high = death
        if high <= q:
            return 'FAIL_DEATH_BY_Q'
        if low <= q:
            return 'UNKNOWN_DEATH_OVERLAPS_Q'
        return 'PASS_Q'
    if admission['destination'] == 'Overleden':
        return 'UNKNOWN_DEATH_TYPE_NO_DATE'
    return 'PASS_Q' if survival_record_after(admission, patient_rows, q) else 'UNKNOWN_Q_SURVIVAL'

def medication_states(rows: list[dict[str, str]], q: int, h: int, c: int) -> tuple[list[int], list[int]]:
    definite: list[int] = []
    unresolved_possible_starts: list[int] = []
    upper = min(h, c)
    for row in rows:
        kind = technical.classify(row)
        if kind.startswith('noncontinuous') or kind == 'definite_nonpositive':
            continue
        start, stop = (integer(row['start']), integer(row['stop']))
        if start is None or stop is None or start >= stop:
            continue
        if start is not None and start > upper:
            continue
        if stop is not None and stop <= q:
            continue
        if start >= q and start <= h and (start < c):
            if kind.startswith('definite_positive_'):
                definite.append(start)
            elif kind.startswith('unresolved'):
                unresolved_possible_starts.append(start)
        elif kind.startswith('unresolved') and start < upper and (stop > q):
            unresolved_possible_starts.append(max(q, start))
    return (sorted(definite), sorted(unresolved_possible_starts))

def restart_binary(rows: list[dict[str, str]], admission: dict[str, str], patient_rows: list[dict[str, str]], t0: int, horizon_hours: int, death: tuple[int, int] | None, conflict: bool) -> tuple[str, int | None, bool]:
    q, h = (t0 + 4 * HOUR, t0 + horizon_hours * HOUR)
    c = integer(admission['dischargedat'])
    if c is None or conflict:
        return ('UNKNOWN_SOURCE_CONFLICT', None, False)
    starts, unresolved_starts = medication_states(rows, q, h, c)
    if death is not None:
        low, high = death
        unresolved_starts = [start for start in unresolved_starts if start < high]
        before_death = [start for start in starts if start < low]
        if before_death:
            return ('YES', before_death[0], any((start <= before_death[0] for start in unresolved_starts)))
        if any((start < high for start in starts)):
            return ('UNKNOWN_DEATH_RESTART_ORDER', None, bool(unresolved_starts))
    elif starts:
        first = starts[0]
        if survival_record_after(admission, patient_rows, first):
            return ('YES', first, any((start <= first for start in unresolved_starts)))
        return ('UNKNOWN_RESTART_SURVIVAL', None, bool(unresolved_starts))
    if unresolved_starts:
        return ('UNKNOWN_MEDICATION_STATE', None, True)
    end_type = destination_type(admission['destination'])
    if c < h:
        if end_type == 'RECORDED_NONDEATH_END':
            return ('NO', None, False)
        if death is not None and death[0] > q and (death[1] <= c):
            return ('NO', None, False)
        return ('UNKNOWN_CARE_END', None, False)
    if death is not None and death[0] > q:
        return ('NO', None, False)
    if c == h and end_type == 'RECORDED_NONDEATH_END':
        return ('NO', None, False)
    if survival_record_after(admission, patient_rows, h):
        return ('NO', None, False)
    return ('UNKNOWN_OBSERVATION_AT_H', None, False)

def death_binary(admission: dict[str, str], patient_rows: list[dict[str, str]], t0: int, horizon_hours: int, death: tuple[int, int] | None, conflict: bool) -> str:
    if conflict:
        return 'UNKNOWN_SOURCE_CONFLICT'
    h = t0 + horizon_hours * HOUR
    if death is not None:
        low, high = death
        if low > t0 and high <= h:
            return 'YES'
        if low > h:
            return 'NO'
        return 'UNKNOWN_DATE_PRECISION'
    if destination_type(admission['destination']) == 'DEATH_TYPE':
        return 'UNKNOWN_DEATH_TYPE_NO_DATE'
    return 'NO' if survival_record_after(admission, patient_rows, h) else 'UNKNOWN_SURVIVAL_AT_H'

def composite_binary(rows: list[dict[str, str]], admission: dict[str, str], patient_rows: list[dict[str, str]], t0: int, horizon_hours: int, death: tuple[int, int] | None, conflict: bool) -> str:
    """Independently adjudicate any definite restart or death by this horizon.

    This is event-source logic, never a sum of the marginal event counts.
    """
    h = t0 + horizon_hours * HOUR
    if conflict:
        return 'UNKNOWN'
    if death is not None and death[0] > t0 and (death[1] <= h):
        return 'YES'
    restart_status, _, _ = restart_binary(rows, admission, patient_rows, t0, horizon_hours, death, conflict)
    if restart_status == 'YES':
        return 'YES'
    death_status = death_binary(admission, patient_rows, t0, horizon_hours, death, conflict)
    if restart_status == 'NO' and death_status == 'NO':
        return 'NO'
    return 'UNKNOWN'

def first_transition(restart: str, first_restart: int | None, admission: dict[str, str], t0: int, death: tuple[int, int] | None, horizon_hours: int) -> str:
    h = t0 + horizon_hours * HOUR
    c = integer(admission['dischargedat'])
    if restart == 'YES' and first_restart is not None:
        return 'RESTART'
    if restart != 'NO' or c is None:
        return 'UNKNOWN_FIRST_TRANSITION'
    end_type = destination_type(admission['destination'])
    if c < h and end_type == 'RECORDED_NONDEATH_END':
        if death is None or death[0] > c:
            return 'NORMAL_UNIT_END'
        if death[1] <= c:
            return 'DEATH_BEFORE_RESTART'
        return 'UNKNOWN_FIRST_TRANSITION'
    if death is not None and death[0] > t0 and (death[1] <= min(c, h)):
        return 'DEATH_BEFORE_RESTART'
    if c >= h and (death is None or death[0] > h):
        return 'STILL_OBSERVED_NO_RESTART'
    return 'UNKNOWN_FIRST_TRANSITION'

def support_status(rows: list[dict[str, str]], itemid: str, t0: int) -> str:
    unresolved = False
    for row in rows:
        if row['itemid'] != itemid:
            continue
        start, stop = (integer(row['start']), integer(row['stop']))
        if start is not None and start > t0:
            continue
        if stop is not None and stop <= t0:
            continue
        if start is None or stop is None or start >= stop:
            unresolved = True
        elif start <= t0 < stop:
            return 'RECORDED_ACTIVE_PROCESS'
    return 'UNKNOWN_PROCESS_INTERVAL' if unresolved else 'NO_MAPPED_PROCESS_RECORDED'

def choose_first_candidate(rows: list[dict[str, str]], admission: dict[str, str], patient_rows: list[dict[str, str]], death: tuple[int, int] | None, conflict: bool, flow: Counter[str]) -> tuple[dict[str, str], int] | None:
    for row in sorted(rows, key=lambda value: (integer(value['t0_ms']), integer(value['episode_index']))):
        proxy = row['preoutcome_proxy_flag']
        if proxy.startswith('FAIL'):
            flow['candidate_proxy_fail'] += 1
            continue
        if proxy.startswith('UNKNOWN'):
            flow['candidate_proxy_unknown_first_block'] += 1
            return None
        if proxy != 'PASS_PREOUTCOME_PROXY_ONLY':
            raise RuntimeError('candidate_proxy_unexpected_state')
        t0 = integer(row['t0_ms'])
        if t0 is None:
            raise RuntimeError('candidate_t0_missing')
        q_status = q_qualification(admission, patient_rows, t0, death, conflict)
        flow[f'candidate_{q_status.lower()}'] += 1
        if q_status.startswith('FAIL'):
            continue
        if q_status.startswith('UNKNOWN'):
            flow['admission_q_unknown_first_block'] += 1
            return None
        return (row, t0)
    return None

def main() -> None:
    check_inputs()
    admissions, patients = read_admissions()
    vaso = read_gzip_by_admission(VASO, 295765)
    process = read_gzip_by_admission(PROCESS, 43379)
    candidates = read_gzip_by_admission(PROXY, 15533)
    flow = Counter()
    outcomes: dict[str, Counter[str]] = defaultdict(Counter)
    selected = 0
    with gzip.open(PART, 'wt', encoding='utf-8', newline='') as saved:
        writer = None
        for admission_id, rows in candidates.items():
            admission = admissions.get(admission_id)
            if admission is None:
                raise RuntimeError('candidate_admission_missing')
            patient_rows = patients[admission['patientid']]
            death, death_conflict = patient_death(patient_rows)
            conflict = source_conflict(admission, patient_rows, death, death_conflict)
            chosen = choose_first_candidate(rows, admission, patient_rows, death, conflict, flow)
            if chosen is None:
                flow['admission_no_definitely_eligible_event'] += 1
                continue
            row, t0 = chosen
            result = {'patientid': admission['patientid'], 'admissionid': admission_id, 'episode_index': row['episode_index'], 't0_ms': str(t0), 'four_drug_nee_peak_lower': row['peak_nee_lower'], 'four_drug_nee_peak_upper': row['peak_nee_upper'], 'positive_hours_segmentwise': row['positive_hours_segmentwise'], 'recorded_ventilation_process_at_t0': support_status(process.get(admission_id, []), '9328', t0), 'recorded_cvvh_process_at_t0': support_status(process.get(admission_id, []), '12465', t0), 'recorded_hemodialysis_process_at_t0': support_status(process.get(admission_id, []), '16363', t0)}
            for horizon in (24, 72):
                restarted, first_at, first_time_unknown = restart_binary(vaso.get(admission_id, []), admission, patient_rows, t0, horizon, death, conflict)
                died = death_binary(admission, patient_rows, t0, horizon, death, conflict)
                composite = composite_binary(vaso.get(admission_id, []), admission, patient_rows, t0, horizon, death, conflict)
                transition = first_transition(restarted, first_at, admission, t0, death, horizon)
                prefix = f'h{horizon}_'
                result[prefix + 'restart'] = restarted
                result[prefix + 'death'] = died
                result[prefix + 'restart_or_death'] = composite
                result[prefix + 'first_transition'] = transition
                result[prefix + 'first_restart_ms'] = '' if first_at is None or first_time_unknown else str(first_at)
                result[prefix + 'restart_time_unknown'] = '1' if first_time_unknown else '0'
                for key in ('restart', 'death', 'restart_or_death', 'first_transition'):
                    outcomes[prefix + key][result[prefix + key]] += 1
            if writer is None:
                writer = csv.DictWriter(saved, fieldnames=list(result))
                writer.writeheader()
            writer.writerow(result)
            selected += 1
            flow['admission_first_definitely_eligible_event'] += 1
    if selected > 15533 or flow['admission_first_definitely_eligible_event'] != selected:
        raise RuntimeError('selected_cohort_count_invariant_failed')
    for horizon in (24, 72):
        for key in ('restart', 'death', 'restart_or_death', 'first_transition'):
            if sum(outcomes[f'h{horizon}_{key}'].values()) != selected:
                raise RuntimeError('outcome_partition_invariant_failed')
    os.replace(PART, OUT)
    AUDIT.write_text(json.dumps({'scope': 'adapted_four_drug_suspected_infection_pure_ic_first_eligible_event', 'source_version': 'AmsterdamUMCdb_v1.0.2', 'no_strict_three_database_claim': True, 'flow_internal_only': dict(flow), 'outcome_partitions_internal_only': {key: dict(value) for key, value in outcomes.items()}, 'no_models_or_cifs_calculated': True}, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'status': 'complete', 'no_models_or_cifs_calculated': True}))
if __name__ == '__main__':
    main()
