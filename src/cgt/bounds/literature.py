"""Finite lower bounds from later sparse-CGT literature.

Reference:
    - C. Shangguan & G. Ge (2016). "New Bounds on the Number of Tests for
      Disjunct Matrices." IEEE Transactions on Information Theory,
      62(12):7518-7521. DOI: 10.1109/TIT.2016.2614726.
      (shangguan_ge_lower_bound, shangguan_ge_constant_weight_lower_bound,
      shangguan_ge_maximum_weight_lower_bound,
      shangguan_ge_weight_aware_lower_bound)
    - H. A. Inan, P. Kairouz & A. Özgür (2020). "Sparse Combinatorial
      Group Testing." IEEE Transactions on Information Theory,
      66(5):2729-2742. DOI: 10.1109/TIT.2019.2953386.
      (inan_kairouz_ozgur_capacity, inan_kairouz_ozgur_lower_bound)

The functions in this module expose the bounds in the same direction as the
rest of :mod:`cgt.bounds`: they return a necessary number of rows/tests.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, isqrt


def _ceil_sqrt(value: int) -> int:
    root = isqrt(value)
    return root if root * root == value else root + 1


def shangguan_ge_lower_bound(n: int, d: int) -> int:
    """Return the Shangguan--Ge (2016) identity-threshold bound.

    Their Corollary 3.4 states

    ``t(d,n) >= min{((15 + sqrt(33))/24) d^2, n}``.

    The comparison with the irrational constant is evaluated exactly: an
    integer ``m`` satisfies ``24 m >= (15 + sqrt(33)) d^2`` iff
    ``(24 m - 15 d^2)^2 >= 33 d^4`` (with a nonnegative left factor).
    """
    if n < 1 or d < 1:
        raise ValueError("n and d must be positive")

    d2 = d * d
    m = (15 * d2 + 23) // 24
    while 24 * m - 15 * d2 < 0 or (24 * m - 15 * d2) ** 2 < 33 * d2 * d2:
        m += 1
    return min(n, m)


def shangguan_ge_constant_weight_lower_bound(
    n: int, d: int, weight: int
) -> int:
    """Return Shangguan--Ge Theorem 1.1 for a constant-weight code.

    The theorem gives ``m >= min(n, (d+1)^2)`` when every column has
    weight exactly ``d+1``.  A return value of zero records that the
    hypothesis is not satisfied, following the convention used by the
    article's optional finite certificates.
    """
    if n < 1 or d < 1 or weight < 1:
        raise ValueError("n, d, and weight must be positive")
    return min(n, (d + 1) ** 2) if weight == d + 1 else 0


def shangguan_ge_maximum_weight_lower_bound(
    n: int, d: int, maximum_weight: int
) -> int:
    """Return the integer form of Shangguan--Ge Corollary 3.5.

    If every column has weight at most ``floor(5d/3)``, their strict
    inequality ``m > d^2+d+1`` (in the nontrivial regime ``n>m``) yields
    ``m >= min(n, d^2+d+2)``.  Zero denotes an inapplicable hypothesis.
    """
    if n < 1 or d < 1 or maximum_weight < 1:
        raise ValueError("n, d, and maximum_weight must be positive")
    if maximum_weight > 5 * d // 3:
        return 0
    return min(n, d * d + d + 2)


def shangguan_ge_weight_aware_lower_bound(
    n: int, d: int, maximum_weight: int
) -> int:
    """Return the strongest applicable weight-aware SG certificate.

    This combines Theorem 1.1 for constant column weight exactly ``d+1``
    with the integer form of Corollary 3.5 for maximum column weight at most
    ``floor(5d/3)``.  Zero denotes that neither hypothesis applies.
    """
    return max(
        shangguan_ge_constant_weight_lower_bound(n, d, maximum_weight),
        shangguan_ge_maximum_weight_lower_bound(n, d, maximum_weight),
    )


def inan_kairouz_ozgur_capacity(m: int, d: int, maximum_weight: int) -> Fraction:
    """Finite private-set capacity proved by Inan--Kairouz--Ozgur.

    Put ``ell = ceil((maximum_weight - 1) / d)``.  For ``ell >= 2``, the
    counting inequality in the proof of their sparse-codeword theorem is

    ``n <= sum_{i=1}^ell C(m,i)
           + C(m,ell+1)/C((ell-1)d+2,ell+1)``.

    For ``maximum_weight <= d`` individual testing is necessary.  The
    special ``ell=1`` theorem is inverted directly in
    :func:`inan_kairouz_ozgur_lower_bound`.
    """
    if m < 0 or d < 1 or maximum_weight < 1:
        raise ValueError("m must be nonnegative and d and maximum_weight positive")
    if maximum_weight <= d:
        return Fraction(m)
    ell = (maximum_weight - 1 + d - 1) // d
    if ell == 1:
        # The exact special-case consequence is n <= m^2/(d(d+1))
        # whenever m<n; representing it as a capacity is convenient for
        # diagnostics, while inversion handles the alternative m=n exactly.
        return Fraction(m * m, d * (d + 1))
    capacity = sum(comb(m, i) for i in range(1, min(ell, m) + 1))
    if ell + 1 <= m:
        denominator = comb((ell - 1) * d + 2, ell + 1)
        capacity += Fraction(comb(m, ell + 1), denominator)
    return Fraction(capacity)


def inan_kairouz_ozgur_lower_bound(n: int, d: int, maximum_weight: int) -> int:
    """Invert the finite IKÖ maximum-column-weight converse."""
    if n < 1 or d < 1 or maximum_weight < 1:
        raise ValueError("n, d, and maximum_weight must be positive")
    if maximum_weight <= d:
        return n

    ell = (maximum_weight - 1 + d - 1) // d
    if ell == 1:
        return min(n, _ceil_sqrt(n * d * (d + 1)))

    def capacity_at_least(m: int) -> bool:
        """Compare with ``n`` without constructing an enormous full sum."""
        total = 0
        for i in range(1, min(ell, m) + 1):
            total += comb(m, i)
            if total >= n:
                return True
        if ell + 1 <= m:
            denominator = comb((ell - 1) * d + 2, ell + 1)
            return total * denominator + comb(m, ell + 1) >= n * denominator
        return total >= n

    lo, hi = 0, n
    while lo < hi:
        mid = (lo + hi) // 2
        if capacity_at_least(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo
