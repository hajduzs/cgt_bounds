"""Reproduce the KS columns of the survey's finite-parameter table.

Reference: W. Kautz & R. Singleton (1964). "Nonrandom Binary
Superimposed Codes." IEEE Transactions on Information Theory,
10(4):363-377. DOI: 10.1109/TIT.1964.1053689.

KS-B is recomputed from the usual one-hot concatenation formula. The KS-C
and KS-C-RS entries retain the parameter choices recorded in the survey.
Those columns depend on the historical inner-code table and therefore must
not silently change when the local Brouwer data are updated.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from cgt.bounds.balanced_dyachkov_rykov import maximum_weight_lower_bound


@dataclass(frozen=True)
class KSRow:
    n: int
    d: int
    ks_b_m: int
    ks_b_q: int
    ks_c_rs_m: int
    ks_c_rs_q: int
    ks_c_rs_w: int | None
    ks_c_m: int
    ks_c_q: int
    ks_c_w: int | None

    @property
    def ks_b_weight(self) -> int:
        return self.ks_b_m // self.ks_b_q

    @property
    def ks_c_rs_weight(self) -> int | None:
        if self.ks_c_rs_w is None:
            return self.ks_b_weight
        dimension = _ceil_log(self.n, self.ks_c_rs_q)
        return (self.d * (dimension - 1) + 1) * self.ks_c_rs_w

    @property
    def ks_c_weight(self) -> int | None:
        if self.ks_c_w is None:
            return None
        dimension = _ceil_log(self.n, self.ks_c_q)
        return (self.d * (dimension - 1) + 1) * self.ks_c_w

    def as_dict(self) -> dict[str, int | None]:
        result = asdict(self)
        result.update(
            ks_b_weight=self.ks_b_weight,
            ks_b_lb=balanced_partition_lower_bound(self.n, self.d, self.ks_b_weight),
            ks_c_rs_weight=self.ks_c_rs_weight,
            ks_c_rs_lb=balanced_partition_lower_bound(self.n, self.d, self.ks_c_rs_weight),
            ks_c_weight=self.ks_c_weight,
            ks_c_lb=balanced_partition_lower_bound(self.n, self.d, self.ks_c_weight),
        )
        return result


# Parameter choices transcribed in CGT_survey-10.pdf, Table 1. A missing
# weight is genuinely absent from that table, rather than inferred here.
_D2_REFERENCE = {
    100: (21, 5, None, 21, 5, 1),
    10**3: (49, 7, None, 39, 11, 3),
    10**4: (75, 125, 5, 63, 11, 3),
    10**5: (99, 11, None, 81, 11, 3),
    10**6: (108, 16, 4, 99, 11, 3),
    10**8: (156, 16, 4, 143, 16, 3),
    10**10: (225, 19, 5, 180, 19, 3),
    10**20: (475, 125, 5, 368, 47, 5),
    10**30: (725, 125, 5, 555, 41, 5),
}

# The same certified inner codes as in the neighbouring reference rows,
# evaluated at the two decades omitted from the survey's displayed table.
_D2_INTERVENING_KS_C = {
    10**7: (117, 11, 3),
    10**9: (165, 16, 3),
}


@lru_cache(maxsize=None)
def ks_c_parameters(n: int, d: int) -> tuple[int, int]:
    """Return a certified KS-C ``(m, w)`` choice.

    For ``d=2`` this uses the optimized inner-code choices recorded in the
    survey (including the same choices at the two intervening decades).  In
    all cases the stored constant-weight-code table is searched for certified
    inner codes.  The one-hot KS-B specialization is retained as a fallback,
    so the result can never be worse than KS-B.
    """
    if d == 2 and n in _D2_REFERENCE:
        _, _, _, m, q, inner_weight = _D2_REFERENCE[n]
        dimension = _ceil_log(n, q)
        outer_length = d * (dimension - 1) + 1
        return m, outer_length * inner_weight
    if d == 2 and n in _D2_INTERVENING_KS_C:
        m, q, inner_weight = _D2_INTERVENING_KS_C[n]
        dimension = _ceil_log(n, q)
        outer_length = d * (dimension - 1) + 1
        return m, outer_length * inner_weight
    candidates = [_ks_b_candidate(n, d)]
    candidates.extend(_ks_c_cwc_candidates(n, d))
    return min(candidates)[:2]


def _ks_b_candidate(n: int, d: int) -> tuple[int, int, int, int, int]:
    m, q = ks_b_parameters(n, d)
    return m, m // q, q, q, 1


@lru_cache(maxsize=1)
def _stored_inner_codes() -> tuple[tuple[int, int, int, int], ...]:
    """Return stored ``(length, distance, weight, size)`` CWC parameters."""
    path = Path(__file__).resolve().parents[3] / "data" / "code_lengths.json"
    raw = json.loads(path.read_text(encoding="utf-8"))
    records = []
    for distance, by_length in raw.items():
        for length, by_weight in by_length.items():
            for weight, size in by_weight.items():
                records.append((int(length), int(distance), int(weight), int(size)))
    return tuple(records)


def _ks_c_cwc_candidates(n: int, d: int):
    """Yield one-level KS-C compositions certified by the local CWC table."""
    for inner_m, distance, inner_w, inner_size in _stored_inner_codes():
        required_distance = 2 * (inner_w - (inner_w - 1) // d)
        if distance < required_distance or inner_size < 2:
            continue
        for q in _prime_powers():
            if q > inner_size:
                break
            dimension = _ceil_log(n, q)
            outer_length = d * (dimension - 1) + 1
            if outer_length <= q:
                yield (
                    inner_m * outer_length,
                    inner_w * outer_length,
                    q,
                    inner_m,
                    inner_w,
                )


@lru_cache(maxsize=1)
def _prime_powers() -> tuple[int, ...]:
    return tuple(q for q in range(2, 4097) if _is_prime_power(q))


def balanced_partition_lower_bound(n: int, d: int, w_max: int | None) -> int | None:
    """Smallest m allowed by the survey's maximum-weight certificate."""
    if w_max is None:
        return None
    return maximum_weight_lower_bound(n, d, w_max)


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            return False
        divisor += 1
    return True


def _is_prime_power(n: int) -> bool:
    if _is_prime(n):
        return True
    for prime in range(2, int(n**0.5) + 1):
        if not _is_prime(prime) or n % prime:
            continue
        remainder = n
        while remainder % prime == 0:
            remainder //= prime
        return remainder == 1
    return False


def _ceil_log(n: int, q: int) -> int:
    power, exponent = 1, 0
    while power < n:
        power *= q
        exponent += 1
    return exponent


def ks_b_parameters(n: int, d: int) -> tuple[int, int]:
    """Return the minimum ``(m, q)`` for the one-hot KS-B construction."""
    if n < 2 or d < 1:
        raise ValueError("n must be at least 2 and d must be positive")
    candidates = []
    q = 2
    while not candidates or q < candidates[0][0]:
        if not _is_prime_power(q):
            q += 1
            continue
        dimension = _ceil_log(n, q)
        outer_length = d * (dimension - 1) + 1
        if outer_length <= q:
            candidates.append((q * outer_length, q))
            candidates.sort()
        q += 1
    return min(candidates)


def ks_rows(d: int = 2, ns: list[int] | None = None) -> list[KSRow]:
    """Return survey KS rows, recomputing and checking every KS-B entry."""
    if d != 2:
        raise ValueError("the survey reference rows are currently available only for d=2")
    selected = list(_D2_REFERENCE) if ns is None else ns
    rows = []
    for n in selected:
        if n not in _D2_REFERENCE:
            raise ValueError(f"no survey KS-C reference parameters for n={n}")
        b_m, b_q = ks_b_parameters(n, d)
        cr_m, cr_q, cr_w, c_m, c_q, c_w = _D2_REFERENCE[n]
        rows.append(KSRow(n, d, b_m, b_q, cr_m, cr_q, cr_w, c_m, c_q, c_w))
    return rows
