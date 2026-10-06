"""Final database-specific algorithm. All clinical outputs are private and must not be committed."""
from __future__ import annotations
import csv
import hashlib
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '00_setup'))
from paths import release_path, private_output
from typing import Any, Iterable, Mapping
ROOT = Path(__file__).resolve().parent
MAPPING_PATH = release_path('amsterdam_antibiotic_mapping')
MAPPING_SHA256 = 'c507a2427420c221b47e25b0ac9351a8679a1051187b6b9e82fe9c7794c37707'
IV_ORDER_CATEGORIES = frozenset({15, 55, 65})
HOUR_MS = 3600000
DAY_MS = 24 * HOUR_MS

@dataclass(frozen=True)
class InfectionProxyDecision:
    status: str
    reasons: tuple[str, ...]

@dataclass(frozen=True)
class _Antibiotic:
    itemid: int
    rank: int
    start: int
    stop: int
    pre_admission_start: bool

@dataclass(frozen=True)
class _UncertainRow:
    start: int | None
    stop: int | None

    def could_be_active_at(self, instant: int) -> bool:
        return (self.start is None or self.start <= instant) and (self.stop is None or self.stop >= instant)

def _int_or_none(value: Any) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        raw = str(value).strip()
        return int(raw) if raw else None
    except (TypeError, ValueError):
        return None

def _decimal_or_none(value: Any) -> Decimal | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        result = Decimal(str(value).strip())
        return result if result.is_finite() else None
    except (InvalidOperation, ValueError):
        return None

def load_pinned_ranks(path: Path=MAPPING_PATH) -> dict[int, int]:
    """Reject a changed or malformed antibiotic dictionary before adjudication."""
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != MAPPING_SHA256:
        raise RuntimeError('antibiotic_mapping_sha256_mismatch')
    with path.open('r', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, strict=True)
        if reader.fieldnames != ['itemid', 'rank']:
            raise RuntimeError('antibiotic_mapping_schema_mismatch')
        ranks: dict[int, int] = {}
        for row in reader:
            itemid, rank = (int(row['itemid']), int(row['rank']))
            if itemid in ranks or rank not in (1, 2, 3, 4):
                raise RuntimeError('antibiotic_mapping_invalid_entry')
            ranks[itemid] = rank
    if len(ranks) != 40:
        raise RuntimeError('antibiotic_mapping_entry_count_mismatch')
    return ranks

def _prophylaxis_status(row: Mapping[str, Any], itemid: int, start: int, admittedat: int, specialty: str | None, cardiac_surgery_specialties: frozenset[str] | None) -> str:
    """Return excluded, eligible or unknown without inferred clinical labels."""
    if itemid == 6919 and admittedat <= start < admittedat + 4 * DAY_MS:
        return 'excluded'
    if itemid == 7208:
        dose = _decimal_or_none(row.get('dose'))
        unit = _int_or_none(row.get('doseunitid'))
        if dose is None or dose <= 0 or unit not in (9, 10):
            return 'unknown'
        if unit == 10 and dose <= Decimal('250') or (unit == 9 and dose <= Decimal('0.25')):
            return 'excluded'
    if itemid == 7064:
        if cardiac_surgery_specialties is None or specialty is None:
            return 'unknown'
        return 'excluded' if specialty in cardiac_surgery_specialties else 'eligible'
    return 'eligible'

def _signature(rows: list[_Antibiotic], instant: int, *, after: bool) -> tuple[int, int] | None:
    active = {row.itemid: row.rank for row in rows if (row.start <= instant and instant < row.stop if after else row.start < instant and instant <= row.stop)}
    if not active:
        return None
    highest = max(active.values())
    return (highest, sum((rank == highest for rank in active.values())))

def _escalates(before: tuple[int, int] | None, after: tuple[int, int] | None) -> bool:
    if after is None:
        return False
    return before is None or after[0] > before[0] or (after[0] == before[0] and after[1] > before[1])

def evaluate_infection_proxy(antibiotic_rows: Iterable[Mapping[str, Any]], *, episode_start: int, t0: int, admittedat: int, specialty: str | None, source_complete: bool, cardiac_surgery_specialties: frozenset[str] | None=None, mapping_path: Path=MAPPING_PATH) -> InfectionProxyDecision:
    """Return PASS/FAIL/UNKNOWN for the strictly pre-t0 antibiotic proxy.

    ``source_complete`` must be true only after the source and restricted extract
    have passed identity, hash, and admission-level completeness checks.  FAIL
    means only that the *proxy was not observed*.  ``specialty`` is never
    heuristically matched: a time-valid cardiac-surgery mapping must be frozen
    separately and passed as ``cardiac_surgery_specialties``.  No such mapping
    exists in v1 of the operational specification.
    """
    ranks = load_pinned_ranks(mapping_path)
    if not all((isinstance(value, int) and (not isinstance(value, bool)) for value in (episode_start, t0, admittedat))):
        return InfectionProxyDecision('UNKNOWN', ('invalid_anchor_time',))
    if not admittedat <= episode_start <= t0:
        return InfectionProxyDecision('UNKNOWN', ('inconsistent_anchor_order',))
    if cardiac_surgery_specialties is not None and (not isinstance(cardiac_surgery_specialties, frozenset)):
        raise TypeError('cardiac_surgery_specialties_must_be_frozenset_or_none')
    lower = episode_start - 72 * HOUR_MS
    upper = min(episode_start + 24 * HOUR_MS, t0)
    definite: list[_Antibiotic] = []
    uncertain: list[_UncertainRow] = []
    admission_ids: set[str] = set()
    blocked_definite_trigger = False
    for row in antibiotic_rows:
        if 'admissionid' in row and str(row['admissionid']).strip():
            admission_ids.add(str(row['admissionid']).strip())
        itemid = _int_or_none(row.get('itemid'))
        if itemid is None:
            uncertain.append(_UncertainRow(None, None))
            continue
        if itemid not in ranks:
            continue
        start = _int_or_none(row.get('start'))
        stop = _int_or_none(row.get('stop'))
        route = _int_or_none(row.get('ordercategoryid'))
        if start is not None and start > t0:
            continue
        if stop is not None and stop <= admittedat:
            continue
        if route is not None and route not in IV_ORDER_CATEGORIES:
            continue
        if start is None or route is None or stop is None or (stop <= start):
            uncertain.append(_UncertainRow(start, stop if stop is not None and stop > start else None))
            continue
        prophylaxis = _prophylaxis_status(row, itemid, start, admittedat, specialty, cardiac_surgery_specialties)
        if prophylaxis == 'excluded':
            continue
        if prophylaxis == 'unknown':
            uncertain.append(_UncertainRow(start, stop))
            continue
        definite.append(_Antibiotic(itemid, ranks[itemid], start, min(stop, t0 + 1), start < admittedat))
    if len(admission_ids) > 1:
        return InfectionProxyDecision('UNKNOWN', ('mixed_admission_source_rows',))
    if not source_complete:
        return InfectionProxyDecision('UNKNOWN', ('source_completeness_unverified',))
    candidate_times = sorted({row.start for row in definite if not row.pre_admission_start and lower <= row.start <= upper})
    for instant in candidate_times:
        before = _signature(definite, instant, after=False)
        after = _signature(definite, instant, after=True)
        if not _escalates(before, after):
            continue
        if any((row.could_be_active_at(instant) for row in uncertain)):
            blocked_definite_trigger = True
            continue
        return InfectionProxyDecision('PASS', ('definite_pre_t0_new_or_escalated_iv_order',))
    if blocked_definite_trigger or any((row.start is None or lower <= row.start <= upper for row in uncertain)):
        return InfectionProxyDecision('UNKNOWN', ('only_potential_trigger_or_active_regimen_uncertain',))
    return InfectionProxyDecision('FAIL', ('no_eligible_trigger_in_complete_mapped_source',))

def infection_flag(rows: list[dict[str, Any]], episode_start: int, t0: int, admitted: int, specialty: str | None) -> tuple[str, str]:
    """Integration interface for a hash-verified, complete one-admission extract.

    The caller must first verify the pinned source manifest and restricted
    antibiotic extraction.  The five-argument integration API cannot itself
    prove source completeness, so it must not be called on a truncated sample.
    No specialty-to-cardiac-surgery codebook was signed in v1; a lone
    vancomycin candidate therefore remains UNKNOWN.
    """
    decision = evaluate_infection_proxy(rows, episode_start=episode_start, t0=t0, admittedat=admitted, specialty=specialty, source_complete=True)
    return (decision.status, ';'.join(decision.reasons))
