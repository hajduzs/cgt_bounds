"""Finite certificate extracted from Ruszinkó's compression proof.

Reference: M. Ruszinkó (1994). "On the Upper Bound of the Size of the
r-Cover-Free Families." Journal of Combinatorial Theory, Series A,
66(2):302-310. DOI: 10.1016/0097-3165(94)90067-1.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import ceil, comb, floor


@lru_cache(maxsize=None)
def ruszinko_capacity(m: int, d: int) -> Fraction:
    """Upper-bound the family size by the finite compression argument.

    Ruszinko's set-compression procedure removes at most ``floor(d/2)``
    columns and leaves a 2-cover-free family whose column weights are at
    most ``floor(2m/d)``.  Lemma 3.2 of the paper, specialized to two
    cover sets and summed over the remaining weight layers, gives

        n <= floor(d/2) + sum_{v <= floor(2m/d)}
             2 binom(m,ceil(v/2)) / binom(v,ceil(v/2)).

    The deliberately unrounded rational value is a valid (though often
    loose) finite upper bound on the number of columns.
    """
    if m < 0 or d < 2:
        raise ValueError("m must be nonnegative and d must be at least 2")
    cutoff = floor(Fraction(2 * m, d))
    capacity = Fraction(d // 2)
    for weight in range(1, cutoff + 1):
        level = ceil(Fraction(weight, 2))
        capacity += Fraction(
            2 * comb(m, level), comb(weight, level)
        )
    return capacity


@lru_cache(maxsize=None)
def ruszinko_lower_bound(n: int, d: int) -> int:
    """Invert :func:`ruszinko_capacity` in the number of rows."""
    if n < 1 or d < 2:
        raise ValueError("n must be positive and d must be at least 2")
    m = 0
    while ruszinko_capacity(m, d) < n:
        m += 1
    return m


@lru_cache(maxsize=None)
def ruszinko_maximum_weight_capacity(
    m: int, d: int, maximum_weight: int
) -> Fraction:
    """Finite weight-layer capacity from Ruszinko's Lemma 3.2.

    For a layer of weight ``v``, the lemma gives

        d binom(m,ceil(v/d)) / binom(v,ceil(v/d)).

    Summing the disjoint weight layers up to ``maximum_weight`` gives a
    valid upper bound for families with maximum column weight at most that
    value.
    """
    if m < 0 or d < 2 or maximum_weight < 1:
        raise ValueError(
            "m must be nonnegative, d at least 2, and maximum_weight positive"
        )
    capacity = Fraction(0)
    for weight in range(1, min(m, maximum_weight) + 1):
        level = ceil(Fraction(weight, d))
        capacity += Fraction(
            d * comb(m, level), comb(weight, level)
        )
    return capacity


@lru_cache(maxsize=None)
def ruszinko_maximum_weight_lower_bound(
    n: int, d: int, maximum_weight: int
) -> int:
    """Invert the finite weight-dependent Ruszinko layer sum."""
    if n < 1 or d < 2 or maximum_weight < 1:
        raise ValueError(
            "n and maximum_weight must be positive and d at least 2"
        )
    m = 1
    while ruszinko_maximum_weight_capacity(m, d, maximum_weight) < n:
        m += 1
    return m
