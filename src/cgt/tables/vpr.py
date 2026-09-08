"""Porat-Rothschild / Gilbert-Varshamov table parameters.

Reference: E. Porat & A. Rothschild (2011). "Explicit Nonadaptive
Combinatorial Group Testing Schemes." IEEE Transactions on Information
Theory, 57(12):7982-7989. DOI: 10.1109/TIT.2011.2163296.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import ceil, log

from cgt.tables.ks import balanced_partition_lower_bound


@dataclass(frozen=True)
class VPRRow:
    n: int
    d: int
    q: int
    outer_length: int
    m: int
    lower_bound: int

    @property
    def w(self) -> int:
        # One nonzero bit per q-ary outer coordinate after one-hot expansion.
        return self.outer_length

    @property
    def additive_gap(self) -> int:
        return self.m - self.lower_bound

    @property
    def gap_ratio(self) -> float:
        return self.m / self.lower_bound

    def as_dict(self):
        result = asdict(self)
        result.update(w=self.w, additive_gap=self.additive_gap, gap_ratio=self.gap_ratio)
        return result


def qary_entropy(q: int, delta: float) -> float:
    if not 0 < delta < 1:
        raise ValueError("delta must lie strictly between zero and one")
    return (
        delta * log(q - 1, q)
        - delta * log(delta, q)
        - (1 - delta) * log(1 - delta, q)
    )


def vpr_parameters(n: int, d: int, q: int) -> tuple[int, int]:
    """Return ``(m, w)`` from the Porat-Rothschild GV parameter formula."""
    if n < 2 or d < 2 or q < 2:
        raise ValueError("require n >= 2, d >= 2, and q >= 2")
    delta = (d - 1) / d
    if delta >= 1 - 1 / q:
        raise ValueError("the target relative distance is outside the q-ary GV range")
    rate = 1 - qary_entropy(q, delta)
    if rate <= 0:
        raise ValueError("the selected q gives no positive GV rate")
    outer_length = ceil(log(n, q) / rate)
    return q * outer_length, outer_length


# The first row uses q=5 in the source table; all later d=2 rows use q=7.
_D2_SURVEY_Q = {
    100: 5, 10**3: 7, 10**4: 7, 10**5: 7, 10**6: 7,
    10**8: 7, 10**10: 7, 10**20: 7, 10**30: 7,
}

# Same pattern as d=2: the first row uses a smaller q (7), every later d=3
# row uses q=11. docs/tables.tex's d=3 table skips n=1000 (only these 7
# rows are tabulated for d=3), unlike d=2's 9 rows.
_D3_SURVEY_Q = {
    100: 7, 10**4: 11, 10**6: 11, 10**8: 11,
    10**10: 11, 10**20: 11, 10**30: 11,
}

_SURVEY_Q = {2: _D2_SURVEY_Q, 3: _D3_SURVEY_Q}


def vpr_rows(d: int = 2, ns: list[int] | None = None) -> list[VPRRow]:
    if d not in _SURVEY_Q:
        raise ValueError("survey VPR parameter choices are currently available only for d=2 and d=3")
    survey_q = _SURVEY_Q[d]
    selected = list(survey_q) if ns is None else ns
    rows = []
    for n in selected:
        if n not in survey_q:
            raise ValueError(f"no survey VPR parameter choice for n={n}, d={d}")
        q = survey_q[n]
        m, weight = vpr_parameters(n, d, q)
        lower_bound = balanced_partition_lower_bound(n, d, weight)
        rows.append(VPRRow(n, d, q, weight, m, lower_bound))
    return rows
