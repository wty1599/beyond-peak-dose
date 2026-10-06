"""Final database-specific algorithm. All clinical outputs are private and must not be committed."""
from __future__ import annotations
import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable
HOUR_MS = 3600000
PRIORITY = (20081, 6638, 6637, 10442)
WEIGHT_IDS = set(PRIORITY)
LACTATE_IDS = {10053, 6837, 9580}

def finite_number(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None

def integer(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

def weight_interval(group: str) -> tuple[float, float] | None:
    if group == '59-':
        return (0.0, 60.0)
    if group == '110+':
        return (110.0, math.inf)
    try:
        left, right = group.split('-', 1)
        low, high = (int(left), int(right) + 1)
    except (ValueError, TypeError):
        return None
    if low in (60, 70, 80, 90, 100) and high == low + 10:
        return (float(low), float(high))
    return None

@dataclass(frozen=True)
class WeightResult:
    lower_kg: float
    upper_kg: float
    source: str
    conflict: bool

def select_weight(rows: Iterable[dict[str, str]], t0: int, group: str) -> WeightResult | None:
    """Highest-priority latest pre-t0 kg, or interval group on no valid point/conflict."""
    by_source: dict[int, list[tuple[int, float]]] = defaultdict(list)
    for row in rows:
        item = integer(row.get('itemid', ''))
        if item not in WEIGHT_IDS or row.get('unitid') != '12' or row.get('unit', '').casefold() != 'kg':
            continue
        measured, registered = (integer(row.get('measuredat', '')), integer(row.get('registeredat', '')))
        value = finite_number(row.get('value', ''))
        if measured is None or registered is None or value is None or (measured > t0) or (registered > t0) or (not 30.0 <= value <= 300.0):
            continue
        by_source[item].append((measured, value))
    latest_by_source: dict[int, tuple[int, float] | None] = {}
    for source in PRIORITY:
        values = by_source.get(source, [])
        if not values:
            continue
        last_time = max((time for time, _ in values))
        tied = {value for time, value in values if time == last_time}
        if len(tied) != 1:
            latest_by_source[source] = None
        else:
            latest_by_source[source] = (last_time, next(iter(tied)))
    chosen: tuple[int, float] | None = None
    chosen_source: int | None = None
    for source in PRIORITY:
        if source in latest_by_source:
            chosen_source = source
            chosen = latest_by_source[source]
            break
    conflict = chosen_source is not None and chosen is None
    if chosen is not None:
        chosen_value = chosen[1]
        for source, candidate in latest_by_source.items():
            if source == chosen_source:
                continue
            if candidate is None:
                conflict = True
                continue
            other = candidate[1]
            if abs(chosen_value - other) > 0.2 * min(chosen_value, other):
                conflict = True
        if not conflict:
            return WeightResult(chosen_value, chosen_value, f'numeric_{chosen_source}', False)
    grouped = weight_interval(group)
    if grouped is not None:
        return WeightResult(grouped[0], grouped[1], 'weightgroup_interval', conflict)
    return None

def lactate_flag(rows: Iterable[dict[str, str]], episode_start: int, t0: int) -> str:
    low = episode_start - 6 * HOUR_MS
    high = min(episode_start + 6 * HOUR_MS, t0)
    valid_values: list[float] = []
    unresolved_window = False
    for row in rows:
        item = integer(row.get('itemid', ''))
        if item not in LACTATE_IDS:
            continue
        measured, registered = (integer(row.get('measuredat', '')), integer(row.get('registeredat', '')))
        if measured is None:
            if registered is not None and registered <= t0:
                unresolved_window = True
            continue
        if measured < low or measured > high:
            continue
        if registered is None:
            unresolved_window = True
            continue
        if registered > t0:
            continue
        value = finite_number(row.get('value', ''))
        if row.get('unitid') != '97' or row.get('unit', '').casefold() != 'mmol/l' or value is None or (value < 0):
            unresolved_window = True
            continue
        valid_values.append(value)
    if any((value > 2.0 for value in valid_values)):
        return 'PASS'
    if unresolved_window or not valid_values:
        return 'UNKNOWN'
    return 'FAIL'

def rate_component(row: dict[str, str]) -> tuple[float, float] | None:
    """Return (already per-kg, absolute microgram/min) NEE component.

    State adjudication is separate. This only converts a positive-dose row.
    """
    drug = integer(row.get('itemid', ''))
    dose = finite_number(row.get('dose', ''))
    if dose is None or dose <= 0:
        return None
    unit = row.get('doseunitid')
    time = row.get('doserateunitid')
    perkg = row.get('doserateperkg')
    if drug == 7229:
        if (unit, time, perkg) == ('10', '5', '0'):
            return (0.0, dose * 1000.0 / 60.0)
        if (unit, time, perkg) == ('10', '4', '0'):
            return (0.0, dose * 1000.0)
        if (unit, time, perkg) == ('11', '4', '0'):
            return (0.0, dose)
        if (unit, time, perkg) == ('11', '4', '1'):
            return (dose, 0.0)
    if drug == 6818:
        if (unit, time, perkg) == ('10', '5', '0'):
            return (0.0, dose * 1000.0 / 60.0)
        if (unit, time, perkg) == ('11', '5', '0'):
            return (0.0, dose / 60.0)
    if drug == 7179:
        if (unit, time, perkg) == ('10', '5', '0'):
            return (0.0, dose * 1000.0 / 60.0 / 100.0)
        if (unit, time, perkg) == ('11', '4', '1'):
            return (dose / 100.0, 0.0)
    if drug == 19929 and (unit, time, perkg) == ('10', '5', '0'):
        return (0.0, dose * 1000.0 / 60.0 / 10.0)
    return None

def _component_at(component: tuple[float, float], weight: float) -> float:
    fixed, absolute = component
    if absolute == 0:
        return fixed
    if weight == math.inf:
        return fixed
    if weight == 0:
        return math.inf
    return fixed + absolute / weight

def peak_nee_bounds(rows: Iterable[dict[str, str]], episode_start: int, t0: int, weight: WeightResult | None) -> tuple[float, float]:
    """Bound peak at all source interval endpoints, max-vs-sum per drug."""
    from build_technical_episodes import classify
    positive: list[tuple[int, int, str, tuple[float, float] | None]] = []
    boundaries = {episode_start, t0}
    for row in rows:
        if not classify(row).startswith('definite_positive_'):
            continue
        start, stop = (integer(row.get('start', '')), integer(row.get('stop', '')))
        if start is None or stop is None or stop <= episode_start or (start >= t0):
            continue
        start, stop = (max(start, episode_start), min(stop, t0))
        if start >= stop:
            continue
        positive.append((start, stop, row['itemid'], rate_component(row)))
        boundaries.add(start)
        boundaries.add(stop)
    if not positive:
        return (0.0, math.inf)
    if weight is None:
        weight = WeightResult(0.0, math.inf, 'missing', False)
    points = sorted(boundaries)
    peak_lower = 0.0
    peak_upper = 0.0
    for left, right in zip(points, points[1:]):
        if left >= right:
            continue
        by_drug: dict[str, list[tuple[float, float] | None]] = defaultdict(list)
        for start, stop, drug, component in positive:
            if start <= left < stop:
                by_drug[drug].append(component)
        lower_now = 0.0
        upper_now = 0.0
        for components in by_drug.values():
            known = [component for component in components if component is not None]
            if known:
                lower_now += max((_component_at(component, weight.upper_kg) for component in known))
                upper_now += sum((_component_at(component, weight.lower_kg) for component in known))
            if len(known) != len(components):
                upper_now = math.inf
        peak_lower = max(peak_lower, lower_now)
        peak_upper = max(peak_upper, upper_now)
    return (peak_lower, peak_upper)

def threshold_flag(lower: float, upper: float, threshold: float=0.1) -> str:
    if lower >= threshold:
        return 'PASS'
    if upper < threshold:
        return 'FAIL'
    return 'UNKNOWN'
