"""Final database-specific algorithm. All clinical outputs are private and must not be committed."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from collections import Counter
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output
EXPECTED_BYTES = 818558125
EXPECTED_ROWS = 4907269
EXPECTED_SHA256 = '3bf8bdbabf3e7b67d7ad6d4b9389f2fd689d70d28c6052cb8013ad4a3cddcaa9'
EXPECTED_HEADER = 'admissionid,orderid,ordercategoryid,ordercategory,itemid,item,isadditive,isconditional,rate,rateunit,rateunitid,ratetimeunitid,doserateperkg,dose,doseunit,doserateunit,doseunitid,doserateunitid,administered,administeredunit,administeredunitid,action,start,stop,duration,solutionitemid,solutionitem,solutionadministered,solutionadministeredunit,fluidin,iscontinuous'
SELECTED_IDS = {7229, 6818, 7179, 19929}
VASOPRESSIN_TERMS = ('vasopressin', 'vasopressine', 'argipressin', 'argipressine', 'pitressin', 'pitressine', 'vasostrict', 'empressine', 'antidiuret')
ANALOG_TERMS = ('terlipressin', 'terlipressine', 'glypressin')
OUTPUT_COLUMNS = tuple(EXPECTED_HEADER.split(','))

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=release_path('amsterdam_raw') / 'drugitems.csv')
    parser.add_argument('--output-dir', type=Path, default=private_output('amsterdam_cohort'))
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    out = args.output_dir.resolve(strict=True)
    expected_source = (release_path('amsterdam_raw') / 'drugitems.csv').resolve(strict=True)
    expected_out = private_output('amsterdam_cohort')
    if args.source.is_symlink() or args.output_dir.is_symlink():
        raise SystemExit('symlink_path_rejected')
    if source != expected_source or out != expected_out.resolve(strict=True):
        raise SystemExit('unexpected_source_or_output_path')
    if not source.is_file() or source.stat().st_size != EXPECTED_BYTES:
        raise SystemExit('source_identity_size_mismatch')
    digest = hashlib.sha256()
    with source.open('rb') as identity_stream:
        for chunk in iter(lambda: identity_stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    if digest.hexdigest() != EXPECTED_SHA256:
        raise SystemExit('source_sha256_mismatch')
    source_out = out / 'vaso_source.csv.gz'
    alias_out = out / 'vasopressin_alias_source.csv.gz'
    audit_out = out / 'vaso_source_extract_audit.json'
    for target in (source_out, alias_out, audit_out):
        if target.exists():
            raise SystemExit(f'refuse_overwrite:{target.name}')
    csv.field_size_limit(16 * 1024 * 1024)
    row_count = 0
    selected_count: Counter[str] = Counter()
    action_count: Counter[str] = Counter()
    alias_count = 0
    analog_count = 0
    source_tmp = out / 'vaso_source.csv.gz.incomplete'
    alias_tmp = out / 'vasopressin_alias_source.csv.gz.incomplete'
    if source_tmp.exists() or alias_tmp.exists():
        raise SystemExit('previous_incomplete_extract_present')
    try:
        with source.open('r', encoding='cp1252', newline='') as stream, gzip.open(source_tmp, 'wt', encoding='utf-8', newline='') as saved, gzip.open(alias_tmp, 'wt', encoding='utf-8', newline='') as aliases:
            reader = csv.DictReader(stream, strict=True)
            if ','.join(reader.fieldnames or ()) != EXPECTED_HEADER:
                raise RuntimeError('source_header_mismatch')
            writer = csv.DictWriter(saved, fieldnames=OUTPUT_COLUMNS)
            alias_writer = csv.DictWriter(aliases, fieldnames=OUTPUT_COLUMNS)
            writer.writeheader()
            alias_writer.writeheader()
            for row in reader:
                row_count += 1
                if row_count > EXPECTED_ROWS:
                    raise RuntimeError('row_count_exceeded')
                if None in row or any((value is None for value in row.values())):
                    raise RuntimeError('nonrectangular_csv_row')
                try:
                    item_id = int(row['itemid'])
                except (TypeError, ValueError) as exc:
                    raise RuntimeError('invalid_itemid') from exc
                if item_id in SELECTED_IDS:
                    writer.writerow({key: row[key] for key in OUTPUT_COLUMNS})
                    selected_count[str(item_id)] += 1
                    action_count[row['action']] += 1
                elif (item := row['item'].casefold()):
                    if any((term in item for term in ANALOG_TERMS)):
                        analog_count += 1
                    elif any((term in item for term in VASOPRESSIN_TERMS)) or re.search('\\badh\\b', item):
                        alias_writer.writerow({key: row[key] for key in OUTPUT_COLUMNS})
                        alias_count += 1
            if row_count != EXPECTED_ROWS:
                raise RuntimeError('row_count_mismatch')
        os.replace(source_tmp, source_out)
        os.replace(alias_tmp, alias_out)
        audit = {'purpose': 'technical_source_extract_no_clinical_outcomes', 'source': str(source), 'source_bytes': EXPECTED_BYTES, 'source_rows': row_count, 'selected_itemid_counts': dict(sorted(selected_count.items())), 'alias_rows_internal_only': alias_count, 'excluded_analog_rows_internal_only': analog_count, 'selected_action_counts_internal_only': dict(sorted(action_count.items())), 'restricted_outputs': [source_out.name, alias_out.name], 'interpretation': 'Rows are not yet adjudicated as positive infusions, NEE, t0, or restarts.'}
        audit_out.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
        print(json.dumps({'status': 'complete', 'source_rows': row_count, 'restricted_selected_rows': sum(selected_count.values()), 'no_outcomes_calculated': True}))
    finally:
        for tmp in (source_tmp, alias_tmp):
            if tmp.exists():
                tmp.unlink()
if __name__ == '__main__':
    main()
