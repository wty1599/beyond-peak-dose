"""Bounded post-restart *charted* four-drug infusion-time summary.

Result-known exploratory addendum. No row-level output.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
HERE=private_output('amsterdam_burden');ROOT=Path(__file__).resolve().parent
STAGE=private_output('amsterdam_cohort')
COHORT=release_path('amsterdam_cohort_input');VASO=STAGE/'vaso_source.csv.gz'
ADMISSIONS=release_path('amsterdam_admissions')
OUT=HERE/'burden_summary.json';HOUR=3600000;DAY=86400000
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def integer(value: str) -> int | None:
    try:
        return int(value) if value != "" else None
    except ValueError:
        return None


def quantiles(values: list[float]) -> dict[str, float] | None:
    if not values:
        return None
    x = sorted(values)

    def q(p: float) -> float:
        z = (len(x) - 1) * p
        l, u = math.floor(z), math.ceil(z)
        return x[l] + (z - l) * (x[u] - x[l])

    return {"n": len(x), "min": x[0], "q1": q(.25), "median": q(.5),
            "q3": q(.75), "max": x[-1]}


def merged_duration(intervals: list[tuple[int, int]]) -> float:
    intervals.sort()
    merged: list[list[int]] = []
    for start, stop in intervals:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(stop, merged[-1][1])
        else:
            merged.append([start, stop])
    return sum((stop - start) / HOUR for start, stop in merged)


def main() -> None:
    if OUT.exists():
        raise RuntimeError("refuse_overwrite")
    import build_technical_episodes as technical
    cohort: dict[str, dict[str, str]] = {}
    yes_total = 0
    unknown_time = 0
    with gzip.open(COHORT, "rt", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, strict=True):
            if row["h72_restart"] != "YES":
                continue
            yes_total += 1
            if not row["h72_first_restart_ms"]:
                unknown_time += 1
                continue
            cohort[row["admissionid"]] = row
    if yes_total != 206 or len(cohort) + unknown_time != 206:
        raise RuntimeError("restart_count_mismatch")
    admissions: dict[str, dict[str, str]] = {}
    patient_deaths: dict[str, set[str]] = defaultdict(set)
    with ADMISSIONS.open("r", encoding="cp1252", newline="") as f:
        for row in csv.DictReader(f, strict=True):
            if row["admissionid"] in cohort:
                admissions[row["admissionid"]] = row
            if row["dateofdeath"]:
                patient_deaths[row["patientid"]].add(row["dateofdeath"])
    if len(admissions) != len(cohort):
        raise RuntimeError("selected_admission_missing")
    source: dict[str, list[dict[str, str]]] = defaultdict(list)
    count_source = 0
    with gzip.open(VASO, "rt", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, strict=True):
            count_source += 1
            if count_source > 295_765:
                raise RuntimeError("vaso_row_bound_exceeded")
            if row["admissionid"] in cohort:
                source[row["admissionid"]].append(row)
    if count_source != 295_765:
        raise RuntimeError("vaso_row_count_mismatch")

    status = Counter()
    charted_hours: list[float] = []
    available_hours: list[float] = []
    for admission_id, row in cohort.items():
        admission = admissions[admission_id]
        start_window = int(row["h72_first_restart_ms"])
        h72 = int(row["t0_ms"]) + 72 * HOUR
        care_end = integer(admission["dischargedat"])
        if care_end is None:
            raise RuntimeError("selected_care_end_invalid")
        end_window = min(h72, care_end)
        if not (start_window < end_window and start_window >= int(row["t0_ms"]) + 4 * HOUR):
            raise RuntimeError("invalid_observation_window")
        available_hours.append((end_window - start_window) / HOUR)
        if care_end < h72:
            status["boundary_early_care_end"] += 1
            if admission["destination"] == "Overleden":
                status["boundary_death_type_care_end"] += 1
        else:
            status["boundary_h72"] += 1
        if admission["destination"] == "Overleden":
            status["death_type_destination_any"] += 1
        dates = patient_deaths[admission["patientid"]]
        if len(dates) > 1:
            status["patient_death_date_conflict"] += 1
        elif dates:
            d = integer(next(iter(dates)))
            if d is None:
                status["patient_death_date_unparseable"] += 1
            elif d - DAY < end_window and d + DAY > start_window:
                status["possible_death_envelope_intersects_window"] += 1
        intervals: list[tuple[int, int]] = []
        unresolved = False
        definite_first_start = False
        crosses = False
        touches = False
        for drug in source.get(admission_id, []):
            kind = technical.classify(drug)
            a, b = integer(drug["start"]), integer(drug["stop"])
            if a is None or b is None or a >= b or a >= end_window or b <= start_window:
                continue
            if kind.startswith("definite_positive_"):
                if a == start_window:
                    definite_first_start = True
                intervals.append((max(a, start_window), min(b, end_window)))
                if b > end_window:
                    crosses = True
                elif b == end_window:
                    touches = True
            elif kind.startswith("unresolved"):
                unresolved = True
        if not definite_first_start:
            raise RuntimeError("stored_first_restart_not_found_in_source")
        if unresolved:
            status["cumulative_hours_unknown_unresolved_state"] += 1
            status["last_interval_state_unknown"] += 1
            continue
        hours = merged_duration(intervals)
        if not (0 < hours <= (end_window - start_window) / HOUR + 1e-9):
            raise RuntimeError("charted_time_out_of_window")
        charted_hours.append(hours)
        if crosses:
            status["positive_interval_crosses_observation_boundary"] += 1
        elif touches:
            status["positive_interval_ends_at_observation_boundary"] += 1
        else:
            status["last_positive_interval_ends_before_boundary"] += 1
    if sum(status[key] for key in
           ("cumulative_hours_unknown_unresolved_state",
            "positive_interval_crosses_observation_boundary",
            "positive_interval_ends_at_observation_boundary",
            "last_positive_interval_ends_before_boundary")) != len(cohort):
        raise RuntimeError("course_partition_failed")
    result = {
        "source_version": "AmsterdamUMCdb v1.0.2",
        "scope": "charted four-drug re-exposure after definite first restart",
        "restart_yes_total": yes_total,
        "first_restart_time_unknown": unknown_time,
        "first_restart_time_known": len(cohort),
        "charted_positive_hours_determinate": quantiles(charted_hours),
        "available_observation_window_hours": quantiles(available_hours),
        "status_counts": dict(status),
        "not_complete_clinical_course_or_postrestart_peak_NEE": True,
        "no_row_level_fields": True,
    }
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": "complete", "restart_yes": yes_total,
                      "first_time_known": len(cohort),
                      "charted_hours_determinate": len(charted_hours)}))


if __name__ == "__main__":
    main()
