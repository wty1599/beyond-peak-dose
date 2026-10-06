"""Final database-specific algorithm. All clinical outputs are private and must not be committed."""
from __future__ import annotations
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
import time
import zipfile
from collections import Counter
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output
ROOT = Path(__file__).resolve().parent
STAGE = private_output('amsterdam_cohort')
REVIEWED = ROOT / 'numeric_reader.py'
INNER = release_path('amsterdam_numeric_zip')
VASO = STAGE / 'vaso_source.csv.gz'
OUT = STAGE / 'numeric_eligibility_source.csv.gz'
PART = STAGE / 'numeric_eligibility_source.csv.gz.incomplete'
AUDIT = STAGE / 'numeric_eligibility_extract_audit.json'
EXPECTED_VASO_ROWS = 295765
EXPECTED_NUMERIC_ROWS = 977625612
TARGET_IDS = {20081, 6638, 6637, 10442, 10053, 6837, 9580}
PROGRESS_ROWS = 20000000

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def load_reviewed_module():
    spec = importlib.util.spec_from_file_location('reviewed_numeric_reader', REVIEWED)
    if spec is None or spec.loader is None:
        raise RuntimeError('reviewed_reader_missing')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def eligible_admissions() -> set[str]:
    chosen: set[str] = set()
    rows = 0
    with gzip.open(VASO, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        if not {'admissionid', 'itemid'}.issubset(reader.fieldnames or []):
            raise RuntimeError('vaso_extract_schema_mismatch')
        for row in reader:
            rows += 1
            if rows > EXPECTED_VASO_ROWS:
                raise RuntimeError('vaso_row_count_exceeded')
            chosen.add(row['admissionid'])
    if rows != EXPECTED_VASO_ROWS:
        raise RuntimeError('vaso_row_count_mismatch')
    return chosen

def main() -> None:
    for item in (INNER, VASO):
        if item.is_symlink() or not item.is_file():
            raise RuntimeError('restricted_input_missing_or_symlink')
    for target in (OUT, PART, AUDIT):
        if target.exists():
            raise RuntimeError(f'refuse_overwrite:{target.name}')
    reviewed = load_reviewed_module()
    if INNER.stat().st_size != reviewed.STAGED_SIZE or sha256_file(INNER) != reviewed.STAGED_SHA256:
        raise RuntimeError('inner_zip_identity_mismatch')
    admitted = eligible_admissions()
    rows_total = 0
    target_counts: Counter[str] = Counter()
    saved_rows = 0
    started = time.monotonic()
    try:
        with zipfile.ZipFile(INNER, 'r', allowZip64=True) as archive:
            members = archive.infolist()
            if len(members) != 1:
                raise RuntimeError('inner_member_count_mismatch')
            member = members[0]
            if member.filename != reviewed.MEMBER or member.compress_type != zipfile.ZIP_DEFLATED or member.flag_bits & 1 or (member.file_size != reviewed.CSV_BYTES) or (member.compress_size != reviewed.CSV_COMPRESSED_BYTES) or (member.CRC != reviewed.CSV_CRC32):
                raise RuntimeError('inner_member_metadata_mismatch')
            with archive.open(member, 'r') as source, gzip.open(PART, 'wt', encoding='utf-8', newline='') as saved:
                bounded = reviewed.BoundedCRCReader(source)
                with io.TextIOWrapper(io.BufferedReader(bounded, buffer_size=256 * 1024), encoding='cp1252', errors='strict', newline='') as stream:
                    reader = csv.reader(stream, delimiter=',', quotechar='"', strict=True)
                    header = next(reader)
                    if header != reviewed.EXPECTED_HEADER:
                        raise RuntimeError('numeric_schema_mismatch')
                    writer = csv.writer(saved)
                    writer.writerow(header)
                    for row in reader:
                        rows_total += 1
                        if rows_total > EXPECTED_NUMERIC_ROWS or len(row) != 15:
                            raise RuntimeError('numeric_row_count_or_width_exceeded')
                        if rows_total % PROGRESS_ROWS == 0:
                            print(json.dumps({'progress_rows': rows_total, 'elapsed_seconds': round(time.monotonic() - started, 1)}), flush=True)
                        if row[0] in admitted and row[1].isdigit() and (int(row[1]) in TARGET_IDS):
                            writer.writerow(row)
                            saved_rows += 1
                            target_counts[row[1]] += 1
                    if rows_total != EXPECTED_NUMERIC_ROWS:
                        raise RuntimeError('numeric_row_count_mismatch')
                    if bounded.total != reviewed.CSV_BYTES or bounded.crc != reviewed.CSV_CRC32:
                        raise RuntimeError('numeric_csv_crc_or_size_mismatch')
        os.replace(PART, OUT)
        AUDIT.write_text(json.dumps({'purpose': 'eligibility_source_only_no_outcomes', 'numeric_source_sha256_verified': True, 'source_rows': rows_total, 'csv_crc32_verified': True, 'restricted_saved_rows': saved_rows, 'target_itemid_counts_internal_only': dict(target_counts), 'restricted_admission_count': len(admitted), 'elapsed_seconds': round(time.monotonic() - started, 1)}, ensure_ascii=False, indent=2), encoding='utf-8')
        print(json.dumps({'status': 'complete', 'source_rows': rows_total, 'csv_crc32_verified': True, 'no_outcomes_calculated': True}))
    except BaseException:
        raise
if __name__ == '__main__':
    main()
