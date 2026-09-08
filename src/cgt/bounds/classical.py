"""Classical finite lower bounds for disjunct matrices.

Reference:
    - W. Kautz & R. Singleton (1964). "Nonrandom Binary Superimposed
      Codes." IEEE Transactions on Information Theory, 10(4):363-377.
      DOI: 10.1109/TIT.1964.1053689. (kautz_singleton_capacity/lower_bound)
    - A. G. D'yachkov & V. V. Rykov (1982). "Bounds on the Length of
      Disjunctive Codes." Problemy Peredachi Informatsii, 18(3):7-13.
      DOI: none exists (predates the DOI system; not indexed in Crossref
      under any translation venue either). (dyachkov_rykov_capacity/
      lower_bound, c_dr_capacity/lower_bound)
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb

from cgt.bounds.balanced_dyachkov_rykov import sperner_lower_bound


def kautz_singleton_capacity(m: int, d: int) -> int:
    """Right-hand side of the Kautz--Singleton counting inequality."""
    if m < 0 or d < 1:
        raise ValueError("m must be nonnegative and d must be positive")
    return sum(comb(m, j) for j in range(d, m + 1))


def kautz_singleton_lower_bound(n: int, d: int) -> int:
    """Invert ``sum_i binom(n,i) <= sum_j binom(m,j)`` in ``m``."""
    if n < 1 or d < 1:
        raise ValueError("n and d must be positive")
    if d > n // 2:
        raise ValueError("the Kautz--Singleton inequality requires d <= n/2")
    required = sum(comb(n, i) for i in range(1, d + 1))
    m = d
    while kautz_singleton_capacity(m, d) < required:
        m += 1
    return m


def c_dr_capacity(m: int, d: int, maximum_weight: int) -> Fraction:
    """Closed relaxation of the classical DR layer capacity.

    This is the quantity denoted ``C_DR`` in the article.  It is retained for
    reproducibility of older measurements; new DR table rows use the exact
    ``U_DR`` quantity returned by :func:`dyachkov_rykov_capacity`.
    """
    if m < 0 or d < 1:
        raise ValueError("m must be nonnegative and d must be positive")
    maximum_weight = min(maximum_weight, m)
    capacity = Fraction(m)
    for weight in range(d + 1, maximum_weight + 1):
        level = (weight + d - 1) // d
        capacity += Fraction(d * d * comb(m, level), comb(d * level, level))
    return capacity


def dyachkov_rykov_capacity(m: int, d: int, maximum_weight: int) -> Fraction:
    """Exact classical DR weight-layer capacity ``U_DR(m,d,W)``."""
    if m < 0 or d < 1:
        raise ValueError("m must be nonnegative and d must be positive")
    maximum_weight = min(maximum_weight, m)
    capacity = Fraction(m)
    for weight in range(d + 1, maximum_weight + 1):
        level = (weight + d - 1) // d
        capacity += Fraction(
            comb(m, level), comb(weight - 1, level - 1)
        )
    return capacity


def maximum_weight_c_dr_lower_bound(n: int, d: int, maximum_weight: int) -> int:
    """Invert the older closed relaxation ``C_DR``."""
    if n < 1 or d < 1 or maximum_weight < 1:
        raise ValueError("n, d, and maximum_weight must be positive")
    m = 1
    while c_dr_capacity(m, d, maximum_weight) < n:
        m += 1
    return m


def maximum_weight_dyachkov_rykov_lower_bound(
    n: int, d: int, maximum_weight: int
) -> int:
    """Invert the exact classical DR layer bound at a supplied weight budget.

    This is the nonrecursive, weight-aware component of the 1982
    D'yachkov--Rykov argument.  Keeping it separate makes it possible to
    compare that classical counting step with the balanced and recursive
    certificates used in the article.
    """
    if n < 1 or d < 1 or maximum_weight < 1:
        raise ValueError("n, d, and maximum_weight must be positive")
    m = 1
    while dyachkov_rykov_capacity(m, d, maximum_weight) < n:
        m += 1
    return m


@lru_cache(maxsize=None)
def dyachkov_rykov_lower_bound(n: int, d: int) -> int:
    """Compute the exact-``U_DR`` classical deletion-recursive bound."""
    if n < 1 or d < 1:
        raise ValueError("n and d must be positive")
    if n == 1:
        return 1
    if d == 1:
        return sperner_lower_bound(n)
    recursive_rows = dyachkov_rykov_lower_bound(n - 1, d - 1)
    m = recursive_rows + 1
    while True:
        maximum_weight = m - recursive_rows
        if dyachkov_rykov_capacity(m, d, maximum_weight) >= n:
            return m
        m += 1


@lru_cache(maxsize=None)
def c_dr_lower_bound(n: int, d: int) -> int:
    """Older deletion recursion based on the closed ``C_DR`` relaxation."""
    if n < 1 or d < 1:
        raise ValueError("n and d must be positive")
    if n == 1:
        return 1
    if d == 1:
        return sperner_lower_bound(n)
    recursive_rows = c_dr_lower_bound(n - 1, d - 1)
    m = recursive_rows + 1
    while c_dr_capacity(m, d, m - recursive_rows) < n:
        m += 1
    return m
