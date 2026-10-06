"""Final database-specific algorithm. All clinical outputs are private and must not be committed."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import os
from collections import Counter
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output
ROOT = Path(__file__).resolve().parent
STAGE = private_output('amsterdam_cohort')
DRUG = release_path('amsterdam_raw') / 'drugitems.csv'
PROCESS = release_path('amsterdam_raw') / 'processitems.csv'
MAP = release_path('amsterdam_antibiotic_mapping')
DRUG_SHA = '3bf8bdbabf3e7b67d7ad6d4b9389f2fd689d70d28c6052cb8013ad4a3cddcaa9'
PROCESS_SHA = 'd40bd9861f9bcf410dc2ae4727715e790b4c426a2007502c370095e91dcec64b'
MAP_SHA = 'c507a2427420c221b47e25b0ac9351a8679a1051187b6b9e82fe9c7794c37707'
DRUG_BYTES = 818558125
PROCESS_BYTES = 13046999
EXPECTED_DRUG_ROWS = 4907269
PROCESS_IDS = {9164, 10430, 20664, 21091, 21092, 9379, 9328, 10740, 12465, 16363, 12634, 12635}

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def checked(path: Path, size: int, sha: str) -> None:
    if path.is_symlink() or not path.is_file() or path.stat().st_size != size:
        raise RuntimeError('source_path_or_size_mismatch')
    if sha256_file(path) != sha:
        raise RuntimeError('source_sha256_mismatch')

def selected_ids() -> set[int]:
    if sha256_file(MAP) != MAP_SHA:
        raise RuntimeError('antibiotic_dictionary_hash_mismatch')
    with MAP.open('r', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        if reader.fieldnames != ['itemid', 'rank']:
            raise RuntimeError('antibiotic_dictionary_schema_mismatch')
        ids = {int(row['itemid']) for row in reader}
    if len(ids) != 40:
        raise RuntimeError('antibiotic_dictionary_count_mismatch')
    return ids

def main() -> None:
    checked(DRUG, DRUG_BYTES, DRUG_SHA)
    checked(PROCESS, PROCESS_BYTES, PROCESS_SHA)
    antibiotics = selected_ids()
    drug_out = STAGE / 'antibiotic_source.csv.gz'
    drug_part = STAGE / 'antibiotic_source.csv.gz.incomplete'
    process_out = STAGE / 'process_source.csv.gz'
    process_part = STAGE / 'process_source.csv.gz.incomplete'
    audit_out = STAGE / 'antibiotic_process_extract_audit.json'
    for target in (drug_out, drug_part, process_out, process_part, audit_out):
        if target.exists():
            raise RuntimeError(f'refuse_overwrite:{target.name}')
    drug_count = 0
    selected_drug = 0
    selected_process = 0
    drug_ids: Counter[str] = Counter()
    process_ids: Counter[str] = Counter()
    with DRUG.open('r', encoding='cp1252', newline='') as source, gzip.open(drug_part, 'wt', encoding='utf-8', newline='') as saved:
        reader = csv.DictReader(source, strict=True)
        if not {'admissionid', 'itemid', 'ordercategoryid', 'start', 'stop', 'dose', 'doseunitid'}.issubset(reader.fieldnames or []):
            raise RuntimeError('drug_schema_mismatch')
        writer = csv.DictWriter(saved, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            drug_count += 1
            if drug_count > EXPECTED_DRUG_ROWS:
                raise RuntimeError('drug_rows_exceeded')
            if row['itemid'].isdigit() and int(row['itemid']) in antibiotics:
                writer.writerow(row)
                selected_drug += 1
                drug_ids[row['itemid']] += 1
    if drug_count != EXPECTED_DRUG_ROWS:
        raise RuntimeError('drug_rows_mismatch')
    os.replace(drug_part, drug_out)
    with PROCESS.open('r', encoding='cp1252', newline='') as source, gzip.open(process_part, 'wt', encoding='utf-8', newline='') as saved:
        reader = csv.DictReader(source, strict=True)
        if reader.fieldnames != ['admissionid', 'itemid', 'item', 'start', 'stop', 'duration']:
            raise RuntimeError('process_schema_mismatch')
        writer = csv.DictWriter(saved, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if row['itemid'].isdigit() and int(row['itemid']) in PROCESS_IDS:
                writer.writerow(row)
                selected_process += 1
                process_ids[row['itemid']] += 1
    os.replace(process_part, process_out)
    audit_out.write_text(json.dumps({'scope': 'outcome_blind_source_extract_only', 'drug_source_sha256_verified': True, 'process_source_sha256_verified': True, 'antibiotic_mapping_sha256_verified': True, 'drug_source_rows': drug_count, 'restricted_antibiotic_rows': selected_drug, 'restricted_process_rows': selected_process, 'drug_item_counts_internal_only': dict(drug_ids), 'process_item_counts_internal_only': dict(process_ids), 'no_outcomes_calculated': True}, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'status': 'complete', 'drug_source_rows': drug_count, 'no_outcomes_calculated': True}))
if __name__ == '__main__':
    main()
