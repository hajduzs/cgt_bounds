"""Füredi's finite bound for r-cover-free families.

Reference: Z. Füredi (1996). "On r-Cover-free Families." Journal of
Combinatorial Theory, Series A, 73(1):172-173. DOI: 10.1006/jcta.1996.0012.
"""

from __future__ import annotations

from math import comb


def furedi_capacity(m: int, d: int) -> int:
    """Return the strict upper threshold in Furedi's theorem.

    Every d-cover-free family on an m-element ground set has cardinality
    strictly smaller than the returned integer.
    """
    if m < 0 or d < 1:
        raise ValueError("m must be nonnegative and d must be positive")
    denominator = comb(d + 1, 2)
    layer = max(0, (m - d + denominator - 1) // denominator)
    return d + comb(m, layer)


def furedi_lower_bound(n: int, d: int) -> int:
    """Minimum row count not excluded by Furedi's strict upper bound."""
    if n < 1 or d < 1:
        raise ValueError("n and d must be positive")
    m = 0
    while n >= furedi_capacity(m, d):
        m += 1
    return m


def furedi_certificate(n: int, d: int) -> dict[str, int | str]:
    m = furedi_lower_bound(n, d)
    denominator = comb(d + 1, 2)
    layer = max(0, (m - d + denominator - 1) // denominator)
    return {
        "name": "Furedi 1996 r-cover-free bound",
        "n": n,
        "d": d,
        "lower_bound": m,
        "binomial_layer": layer,
        "strict_capacity": furedi_capacity(m, d),
    }
