"""Survey rows for the Indyk-Ngo-Rudra (INR) parameter formula.

Reference: P. Indyk, H. Q. Ngo & A. Rudra (2010). "Efficiently
Decodable Non-Adaptive Group Testing." Proc. ACM-SIAM Symposium on
Discrete Algorithms (SODA), 1126-1142. DOI: 10.1137/1.9781611973075.91.

They give a concatenated construction (Reed-Solomon outer code, random
list-disjunct inner code) that is d-isolate and efficiently decodable
with m <= 4800 * d^2 * log2(n) tests. Only the parameter formula is
reproduced here, not the explicit construction.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, log2

from cgt.bounds.furedi import furedi_lower_bound


@dataclass(frozen=True)
class INRRow:
    n: int
    d: int
    m: int
    lower_bound: int

    @property
    def additive_gap(self) -> int:
        return self.m - self.lower_bound

    @property
    def gap_ratio(self) -> float:
        return self.m / self.lower_bound

    def as_dict(self):
        return {
            "n": self.n,
            "d": self.d,
            "m": self.m,
            "lower_bound": self.lower_bound,
            "additive_gap": self.additive_gap,
            "gap_ratio": self.gap_ratio,
        }


SURVEY_N = [100, 10**3, 10**4, 10**5, 10**6, 10**8, 10**10, 10**20, 10**30]


def inr_parameters(n: int, d: int) -> int:
    """Return ``m`` from the INR formula, m <= 4800 * d^2 * log2(n)."""
    if n < 2 or d < 2:
        raise ValueError("require n >= 2 and d >= 2")
    return ceil(4800 * d * d * log2(n))


def inr_rows(d: int = 2, ns: list[int] | None = None) -> list[INRRow]:
    selected = SURVEY_N if ns is None else ns
    rows = []
    for n in selected:
        m = inr_parameters(n, d)
        lower_bound = furedi_lower_bound(n, d)
        rows.append(INRRow(n, d, m, lower_bound))
    return rows
