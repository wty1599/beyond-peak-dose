"""Final database-specific algorithm. All clinical outputs are private and must not be committed."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import os
import math
from collections import Counter, defaultdict
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output
ROOT = Path(__file__).resolve().parent
STAGE = private_output('amsterdam_cohort')
ADMISSIONS = release_path('amsterdam_raw') / 'admissions.csv'
VASO = STAGE / 'vaso_source.csv.gz'
OUT = STAGE / 'technical_episodes_v4.csv.gz'
PART = STAGE / 'technical_episodes_v4.csv.gz.incomplete'
AUDIT = STAGE / 'technical_episodes_v4_audit.json'
EXPECTED_ADMISSIONS_SIZE = 2855342
EXPECTED_ADMISSIONS_SHA = '6ecd5b759a88e3db83e5f06ccc1211549ed345e83c29be38b6623ae6517aae25'
EXPECTED_VASO_ROWS = 295765
FOUR_HOURS_MS = 240 * 60000
VALID_DOSING = {('7229', '10', '5', '0'), ('7229', '10', '4', '0'), ('7229', '11', '4', '0'), ('7229', '11', '4', '1'), ('6818', '10', '5', '0'), ('6818', '11', '5', '0'), ('7179', '10', '5', '0'), ('7179', '11', '4', '1'), ('19929', '10', '5', '0')}
VALID_FLUID_RATE = ('6', '5')
VALID_ADMINISTERED_UNIT_IDS = {'10', '11'}
ADULT_GROUPS = {'18-39', '40-49', '50-59', '60-69', '70-79', '80+'}

def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def integer(value: str) -> int | None:
    if not value or len(value) > 20:
        return None
    try:
        number = int(value)
    except ValueError:
        return None
    return number if abs(number) < 10 ** 16 else None

def positive_number(value: str) -> bool | None:
    if not value:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number > 0

def admissions_by_id() -> dict[str, tuple[str, str, str, int | None, int | None]]:
    if ADMISSIONS.is_symlink() or not ADMISSIONS.is_file() or ADMISSIONS.stat().st_size != EXPECTED_ADMISSIONS_SIZE or (file_sha256(ADMISSIONS) != EXPECTED_ADMISSIONS_SHA):
        raise RuntimeError('admissions_identity_mismatch')
    result = {}
    with ADMISSIONS.open('r', encoding='cp1252', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        needed = {'admissionid', 'patientid', 'location', 'agegroup', 'admittedat', 'dischargedat'}
        if not needed.issubset(reader.fieldnames or []):
            raise RuntimeError('admissions_schema_mismatch')
        for row in reader:
            identifier = row['admissionid']
            if identifier in result:
                raise RuntimeError('duplicate_admission_id')
            result[identifier] = (row['patientid'], row['location'], row['agegroup'], integer(row['admittedat']), integer(row['dischargedat']))
    if len(result) != 23106:
        raise RuntimeError('admissions_row_count_mismatch')
    return result

def classify(row: dict[str, str]) -> str:
    if row['action'].casefold() == 'flush':
        return 'noncontinuous_flush'
    dosing = (row['itemid'], row['doseunitid'], row['doserateunitid'], row['doserateperkg'])
    if row['itemid'] == '6818' and (not row['doserateunitid']) and (row['iscontinuous'] == '0') and (positive_number(row['rate']) is not True):
        return 'noncontinuous_bolus_like'
    rate = positive_number(row['rate'])
    administered = positive_number(row['administered'])
    dose = positive_number(row['dose'])
    if rate is True and administered is True:
        if (row['rateunitid'], row['ratetimeunitid']) != VALID_FLUID_RATE or row['administeredunitid'] not in VALID_ADMINISTERED_UNIT_IDS:
            return 'unresolved_fluid_or_amount_unit'
        return 'definite_positive_nee_known' if dosing in VALID_DOSING and dose is True else 'definite_positive_nee_unknown'
    if rate is False and administered is False and (dose is False):
        return 'definite_nonpositive'
    return 'unresolved_sign'

def episode_rows(intervals: list[tuple[int, int]]) -> list[tuple[int, int, float, int]]:
    """Return (start,t0,segmentwise_positive_hours,number_of_intervals)."""
    intervals.sort()
    grouped: list[list[tuple[int, int]]] = []
    max_end: int | None = None
    for start, end in intervals:
        if max_end is None or start >= max_end + FOUR_HOURS_MS:
            grouped.append([])
            max_end = None
        grouped[-1].append((start, end))
        max_end = end if max_end is None else max(max_end, end)
    result = []
    for group in grouped:
        merged: list[list[int]] = []
        for start, end in group:
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        result.append((group[0][0], max((end for _, end in group)), sum(((end - start) / 3600000.0 for start, end in merged)), len(group)))
    return result

def main() -> None:
    for target in (OUT, PART, AUDIT):
        if target.exists():
            raise RuntimeError(f'refuse_overwrite:{target.name}')
    admissions = admissions_by_id()
    by_admission: dict[str, list[tuple[int, int]]] = defaultdict(list)
    unresolved: dict[str, list[tuple[int, int]]] = defaultdict(list)
    row_classes: Counter[str] = Counter()
    location_classes: Counter[str] = Counter()
    age_classes: Counter[str] = Counter()
    rows = 0
    with gzip.open(VASO, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        for row in reader:
            rows += 1
            if rows > EXPECTED_VASO_ROWS:
                raise RuntimeError('vaso_source_row_count_exceeded')
            admission = admissions.get(row['admissionid'])
            if admission is None:
                row_classes['missing_admission'] += 1
                continue
            _, location, agegroup, admitted, discharged = admission
            location_classes[location] += 1
            if location != 'IC':
                continue
            if agegroup not in ADULT_GROUPS:
                age_classes[agegroup or 'MISSING'] += 1
                continue
            start, stop = (integer(row['start']), integer(row['stop']))
            if admitted is None or discharged is None or admitted >= discharged:
                row_classes['invalid_admission_bounds'] += 1
                continue
            if start is None or stop is None or start >= stop:
                row_classes['invalid_order_interval'] += 1
                continue
            start, stop = (max(start, admitted), min(stop, discharged))
            if start >= stop:
                row_classes['outside_unit_bounds'] += 1
                continue
            category = classify(row)
            row_classes[category] += 1
            if category.startswith('definite_positive_'):
                by_admission[row['admissionid']].append((start, stop))
            elif category.startswith('unresolved'):
                unresolved[row['admissionid']].append((start, stop))
    if rows != EXPECTED_VASO_ROWS:
        raise RuntimeError('vaso_source_row_count_mismatch')
    saved = 0
    with gzip.open(PART, 'wt', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['patientid', 'admissionid', 'episode_index', 'episode_start_ms', 't0_ms', 'positive_hours_segmentwise', 'positive_interval_rows', 'unresolved_overlap_episode_or_q', 'unresolved_prior_split'])
        for admission_id, intervals in sorted(by_admission.items(), key=lambda pair: int(pair[0])):
            patientid = admissions[admission_id][0]
            prior_t0: int | None = None
            for index, (start, t0, positive_hours, n_rows) in enumerate(episode_rows(intervals), 1):
                unresolved_close = any((u_start < t0 + FOUR_HOURS_MS and u_end > start for u_start, u_end in unresolved.get(admission_id, ())))
                unresolved_prior = any((u_start < start and u_end > start - FOUR_HOURS_MS for u_start, u_end in unresolved.get(admission_id, ())))
                if prior_t0 is not None:
                    unresolved_prior = unresolved_prior or any((u_start < start and u_end > prior_t0 for u_start, u_end in unresolved.get(admission_id, ())))
                writer.writerow([patientid, admission_id, index, start, t0, positive_hours, n_rows, int(unresolved_close), int(unresolved_prior)])
                saved += 1
                prior_t0 = t0
    os.replace(PART, OUT)
    AUDIT.write_text(json.dumps({'purpose': 'technical_episodes_only_no_clinical_eligibility_or_outcomes', 'source_rows': rows, 'technical_episodes_internal_only': saved, 'class_counts_internal_only': dict(row_classes), 'location_row_counts_internal_only': dict(location_classes), 'unknown_age_source_rows_internal_only': dict(age_classes), 'no_outcomes_calculated': True}, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'status': 'complete', 'source_rows': rows, 'no_outcomes_calculated': True}))
if __name__ == '__main__':
    main()
