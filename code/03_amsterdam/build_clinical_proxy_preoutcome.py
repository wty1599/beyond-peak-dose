"""Final database-specific algorithm. All clinical outputs are private and must not be committed."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output
import build_technical_episodes as technical
from dose_weight_proxy import HOUR_MS, lactate_flag, peak_nee_bounds, select_weight, threshold_flag
from infection_proxy_preoutcome import infection_flag
ROOT = Path(__file__).resolve().parent
STAGE = private_output('amsterdam_cohort')
ADMISSIONS = release_path('amsterdam_raw') / 'admissions.csv'
TECHNICAL = STAGE / 'technical_episodes_v4.csv.gz'
TECHNICAL_AUDIT = STAGE / 'technical_episodes_v4_audit.json'
VASO = STAGE / 'vaso_source.csv.gz'
NUMERIC = STAGE / 'numeric_eligibility_source.csv.gz'
ANTIBIOTIC = STAGE / 'antibiotic_source.csv.gz'
PROCESS = STAGE / 'process_source.csv.gz'
RAW_SOURCE_SHA = {release_path('amsterdam_raw') / 'drugitems.csv': '3bf8bdbabf3e7b67d7ad6d4b9389f2fd689d70d28c6052cb8013ad4a3cddcaa9', release_path('amsterdam_raw') / 'processitems.csv': 'd40bd9861f9bcf410dc2ae4727715e790b4c426a2007502c370095e91dcec64b', release_path('amsterdam_raw') / 'numericitems.zip': '16bf18e7d4261fc69682bc0c8f800a9b9f986c2389e7f4d3fe13e3036909a3fb'}
EXTRACT_AUDITS = {VASO: (STAGE / 'vaso_source_extract_audit.json', 'selected_itemid_counts', 295765, 'source_rows', 4907269), NUMERIC: (STAGE / 'numeric_eligibility_extract_audit.json', 'restricted_saved_rows', 165319, 'source_rows', 977625612), ANTIBIOTIC: (STAGE / 'antibiotic_process_extract_audit.json', 'restricted_antibiotic_rows', 228728, 'drug_source_rows', 4907269), PROCESS: (STAGE / 'antibiotic_process_extract_audit.json', 'restricted_process_rows', 43379, 'drug_source_rows', 4907269)}
OUT = STAGE / 'clinical_proxy_candidates_v1.csv.gz'
PART = STAGE / 'clinical_proxy_candidates_v1.csv.gz.incomplete'
AUDIT = STAGE / 'clinical_proxy_candidates_v1_audit.json'
ADMISSIONS_BYTES = 2855342
ADMISSIONS_SHA = '6ecd5b759a88e3db83e5f06ccc1211549ed345e83c29be38b6623ae6517aae25'
MCS_IDS = {'9164', '10430', '20664'}
FOUR_HOURS = 4 * HOUR_MS

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def integer(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

def records_by_admission(path: Path, expected: set[str], expected_rows: int) -> tuple[dict[str, list[dict[str, str]]], str]:
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f'restricted_input_missing:{path.name}')
    result: dict[str, list[dict[str, str]]] = defaultdict(list)
    count = 0
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        if not expected.issubset(reader.fieldnames or []):
            raise RuntimeError(f'restricted_input_schema_mismatch:{path.name}')
        for row in reader:
            count += 1
            if count > expected_rows:
                raise RuntimeError(f'restricted_input_rows_exceeded:{path.name}')
            result[row['admissionid']].append(row)
    if count != expected_rows:
        raise RuntimeError(f'restricted_input_rows_mismatch:{path.name}')
    digest = sha256_file(path)
    return (result, digest)

def validate_extract_audit(path: Path) -> int:
    audit_path, count_field, expected_count, source_field, expected_source = EXTRACT_AUDITS[path]
    if audit_path.is_symlink() or not audit_path.is_file():
        raise RuntimeError(f'extract_audit_missing:{path.name}')
    audit = json.loads(audit_path.read_text(encoding='utf-8'))
    if count_field == 'selected_itemid_counts':
        actual_count = sum(audit.get(count_field, {}).values())
    else:
        actual_count = audit.get(count_field)
    if actual_count != expected_count or audit.get(source_field) != expected_source:
        raise RuntimeError(f'extract_audit_count_mismatch:{path.name}')
    if path == NUMERIC and (not (audit.get('numeric_source_sha256_verified') and audit.get('csv_crc32_verified'))):
        raise RuntimeError('numeric_source_not_verified')
    if path in (ANTIBIOTIC, PROCESS) and (not (audit.get('drug_source_sha256_verified') and audit.get('process_source_sha256_verified') and audit.get('antibiotic_mapping_sha256_verified'))):
        raise RuntimeError('antibiotic_process_source_not_verified')
    return expected_count

def load_admissions() -> dict[str, dict[str, str]]:
    if ADMISSIONS.is_symlink() or not ADMISSIONS.is_file() or ADMISSIONS.stat().st_size != ADMISSIONS_BYTES or (sha256_file(ADMISSIONS) != ADMISSIONS_SHA):
        raise RuntimeError('admissions_source_identity_mismatch')
    result = {}
    with ADMISSIONS.open('r', encoding='cp1252', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        needed = {'patientid', 'admissionid', 'location', 'agegroup', 'admittedat', 'weightgroup', 'specialty'}
        if not needed.issubset(reader.fieldnames or []):
            raise RuntimeError('admissions_schema_mismatch')
        for row in reader:
            if row['admissionid'] in result:
                raise RuntimeError('duplicate_admission_id')
            result[row['admissionid']] = {field: row[field] for field in needed}
    if len(result) != 23106:
        raise RuntimeError('admissions_count_mismatch')
    return result

def mcs_flag(rows: list[dict[str, str]], start: int, t0: int) -> str:
    unresolved = False
    for row in rows:
        if row['itemid'] not in MCS_IDS:
            continue
        left, right = (integer(row['start']), integer(row['stop']))
        if left is not None and left >= t0:
            continue
        if right is not None and right <= start:
            continue
        if left is None or right is None or left >= right:
            unresolved = True
            continue
        if left < t0 and right > start:
            return 'FAIL_RECORDED_MCS'
    return 'UNKNOWN_MCS_SOURCE' if unresolved else 'PASS_NO_MAPPED_MCS'

def make_preoutcome_flags(technical_row: dict[str, str], admission: dict[str, str], vaso: list[dict[str, str]], numeric: list[dict[str, str]], antibiotics: list[dict[str, str]], process: list[dict[str, str]]) -> dict[str, str]:
    start, t0 = (integer(technical_row['episode_start_ms']), integer(technical_row['t0_ms']))
    admitted = integer(admission['admittedat'])
    hours = float(technical_row['positive_hours_segmentwise'])
    if start is None or t0 is None or admitted is None or (not math.isfinite(hours)):
        raise RuntimeError('invalid_technical_or_admission_time')
    q = t0 + FOUR_HOURS
    weight = select_weight(numeric, t0, admission['weightgroup'])
    peak_low, peak_high = peak_nee_bounds(vaso, start, t0, weight)
    duration = 'PASS' if hours >= 6.0 else 'FAIL'
    peak = threshold_flag(peak_low, peak_high)
    lactate = lactate_flag(numeric, start, t0)
    infection, infection_reason = infection_flag(antibiotics, start, t0, admitted, admission['specialty'])
    mcs = mcs_flag(process, start, t0)
    unresolved_all = []
    for row in vaso:
        if not technical.classify(row).startswith('unresolved'):
            continue
        left, right = (integer(row['start']), integer(row['stop']))
        if left is not None and right is not None and (left < right):
            unresolved_all.append((left, right))
    uncertain_episode = technical_row['unresolved_overlap_episode_or_q'] == '1' or technical_row['unresolved_prior_split'] == '1'
    unresolved_before_candidate = any((left < start for left, _ in unresolved_all))
    positive_in_q = any((technical.classify(row).startswith('definite_positive_') and (left := integer(row['start'])) is not None and ((right := integer(row['stop'])) is not None) and (left < q) and (right > t0) for row in vaso))
    q_drug = 'FAIL_RECORDED_FOUR_DRUG_BEFORE_Q' if positive_in_q else 'UNKNOWN_UNRESOLVED_DRUG' if uncertain_episode or any((left < q and right > t0 for left, right in unresolved_all)) else 'PASS_NO_RECORDED_FOUR_DRUG_RESTART_BEFORE_Q'
    if admission['location'] != 'IC' or admission['agegroup'] not in technical.ADULT_GROUPS:
        proxy = 'UNKNOWN_ADMISSION_SCOPE'
    elif uncertain_episode or unresolved_before_candidate:
        proxy = 'UNKNOWN_CLINICAL_PROXY_OR_FIRST_IDENTITY'
    elif mcs == 'FAIL_RECORDED_MCS' or duration == 'FAIL' or peak == 'FAIL' or (lactate == 'FAIL') or (infection == 'FAIL') or q_drug.startswith('FAIL'):
        proxy = 'FAIL_CLINICAL_PROXY'
    elif q_drug.startswith('UNKNOWN') or mcs.startswith('UNKNOWN') or peak == 'UNKNOWN' or (lactate == 'UNKNOWN') or (infection == 'UNKNOWN'):
        proxy = 'UNKNOWN_CLINICAL_PROXY_OR_FIRST_IDENTITY'
    else:
        proxy = 'PASS_PREOUTCOME_PROXY_ONLY'
    return {'patientid': technical_row['patientid'], 'admissionid': technical_row['admissionid'], 'episode_index': technical_row['episode_index'], 'episode_start_ms': str(start), 't0_ms': str(t0), 'q_ms': str(q), 'positive_hours_segmentwise': technical_row['positive_hours_segmentwise'], 'duration_flag': duration, 'weight_source': weight.source if weight else 'missing', 'weight_conflict': ('1' if weight.conflict else '0') if weight else 'UNKNOWN_MISSING_OR_CONFLICT', 'peak_nee_lower': str(peak_low), 'peak_nee_upper': str(peak_high), 'peak_flag': peak, 'lactate_flag': lactate, 'infection_flag': infection, 'infection_reason': infection_reason, 'mcs_flag': mcs, 'q_drug_flag': q_drug, 'unresolved_episode': '1' if uncertain_episode else '0', 'unresolved_earlier_state': '1' if unresolved_before_candidate else '0', 'preoutcome_proxy_flag': proxy}

def main() -> None:
    for target in (OUT, PART, AUDIT):
        if target.exists():
            raise RuntimeError(f'refuse_overwrite:{target.name}')
    for raw_path, expected_sha in RAW_SOURCE_SHA.items():
        if raw_path.is_symlink() or not raw_path.is_file() or sha256_file(raw_path) != expected_sha:
            raise RuntimeError(f'raw_source_identity_mismatch:{raw_path.name}')
    if TECHNICAL.is_symlink() or not TECHNICAL.is_file() or TECHNICAL_AUDIT.is_symlink() or (not TECHNICAL_AUDIT.is_file()):
        raise RuntimeError('technical_episode_input_missing')
    technical_audit = json.loads(TECHNICAL_AUDIT.read_text(encoding='utf-8'))
    if technical_audit.get('source_rows') != 295765 or technical_audit.get('no_outcomes_calculated') is not True:
        raise RuntimeError('technical_episode_audit_mismatch')
    admissions = load_admissions()
    source_counts = {path: validate_extract_audit(path) for path in EXTRACT_AUDITS}
    vasopressors, vaso_sha = records_by_admission(VASO, {'admissionid', 'itemid', 'start', 'stop'}, source_counts[VASO])
    numeric, numeric_sha = records_by_admission(NUMERIC, {'admissionid', 'itemid', 'value', 'measuredat', 'registeredat'}, source_counts[NUMERIC])
    antibiotics, antibiotic_sha = records_by_admission(ANTIBIOTIC, {'admissionid', 'itemid', 'start', 'stop', 'ordercategoryid'}, source_counts[ANTIBIOTIC])
    process, process_sha = records_by_admission(PROCESS, {'admissionid', 'itemid', 'start', 'stop'}, source_counts[PROCESS])
    counter: Counter[str] = Counter()
    flags: dict[str, Counter[str]] = defaultdict(Counter)
    rows_total = 0
    with gzip.open(TECHNICAL, 'rt', encoding='utf-8', newline='') as stream, gzip.open(PART, 'wt', encoding='utf-8', newline='') as saved:
        reader = csv.DictReader(stream, strict=True)
        required = {'patientid', 'admissionid', 'episode_index', 'episode_start_ms', 't0_ms', 'positive_hours_segmentwise', 'unresolved_overlap_episode_or_q', 'unresolved_prior_split'}
        if not required.issubset(reader.fieldnames or []):
            raise RuntimeError('technical_episode_v4_schema_mismatch')
        writer = None
        for row in reader:
            rows_total += 1
            if rows_total > 100000:
                raise RuntimeError('technical_episode_row_bound_exceeded')
            admission = admissions.get(row['admissionid'])
            if admission is None:
                raise RuntimeError('technical_episode_missing_admission')
            result = make_preoutcome_flags(row, admission, vasopressors.get(row['admissionid'], []), numeric.get(row['admissionid'], []), antibiotics.get(row['admissionid'], []), process.get(row['admissionid'], []))
            if writer is None:
                writer = csv.DictWriter(saved, fieldnames=list(result))
                writer.writeheader()
            writer.writerow(result)
            counter[result['preoutcome_proxy_flag']] += 1
            for field in ('duration_flag', 'peak_flag', 'lactate_flag', 'infection_flag', 'mcs_flag', 'q_drug_flag'):
                flags[field][result[field]] += 1
    if rows_total != technical_audit.get('technical_episodes_internal_only'):
        raise RuntimeError('technical_episode_count_mismatch')
    os.replace(PART, OUT)
    age_unknown = sum((1 for admission in admissions.values() if admission['location'] == 'IC' and admission['agegroup'] not in technical.ADULT_GROUPS))
    AUDIT.write_text(json.dumps({'scope': 'preoutcome_clinical_proxy_no_death_restart_or_care_type', 'technical_episode_rows': rows_total, 'preoutcome_status_counts_internal_only': dict(counter), 'component_counts_internal_only': {key: dict(value) for key, value in flags.items()}, 'source_extract_sha256': {'vaso': vaso_sha, 'numeric': numeric_sha, 'antibiotic': antibiotic_sha, 'process': process_sha, 'technical': sha256_file(TECHNICAL)}, 'pure_ic_admissions_with_unknown_age_internal_only': age_unknown, 'no_outcomes_calculated': True}, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'status': 'complete', 'no_outcomes_calculated': True}))
if __name__ == '__main__':
    main()
