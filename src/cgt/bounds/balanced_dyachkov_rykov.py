"""Balanced D'yachkov-Rykov lower bound.

Reference: this paper (see the accompanying article for the theorem
statement and proof); builds on the classical own-subset counting
argument of A. G. D'yachkov & V. V. Rykov (1982), "Bounds on the Length
of Disjunctive Codes," Problemy Peredachi Informatsii, 18(3):7-13
(DOI: none exists; predates the DOI system, not indexed in Crossref
under any translation venue either).

All capacity calculations use exact rational arithmetic.  In particular,
the summands are not rounded separately before the final comparison.
"""

from __future__ import annotations

from functools import lru_cache
from fractions import Fraction
from math import comb, exp, floor, lgamma

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def balanced_partition_capacity(m: int, d: int, maximum_weight: int) -> Fraction:
    """Return ``U_bal(m, d, maximum_weight)`` exactly."""
    if m < 0 or d < 1:
        raise ValueError("m must be nonnegative and d must be positive")
    if maximum_weight <= 0:
        return Fraction(0)
    maximum_weight = min(maximum_weight, m)
    level = (maximum_weight + d - 1) // d
    remainder = maximum_weight - d * (level - 1)
    capacity = sum(
        (Fraction(d * comb(m, j), comb(d * j, j)) for j in range(1, level)),
        Fraction(0),
    )
    capacity += Fraction(remainder * comb(m, level), comb(maximum_weight, level))
    return capacity


def maximum_weight_lower_bound(n: int, d: int, maximum_weight: int) -> int:
    """Invert the balanced-partition certificate in the number of rows."""
    if n < 1 or d < 1 or maximum_weight < 1:
        raise ValueError("n, d, and maximum_weight must be positive")
    # For m < maximum_weight the constraint is vacuous rather than
    # impossible: every column automatically has weight at most m.
    m = 1
    while balanced_partition_capacity(m, d, maximum_weight) < n:
        m += 1
    return m


@lru_cache(maxsize=None)
def average_weight_capacity(
    m: int, d: int, average_weight: Fraction | int
) -> Fraction:
    """Return the average-weight capacity from the article's truncation bound.

    The direct ``W=m`` maximum-weight certificate is included separately.
    This is both sharper than applying the Markov factor at ``W >= m`` and
    handles the limiting case in the displayed infimum over truncation levels.
    """
    average_weight = Fraction(average_weight)
    if m < 0 or d < 1 or average_weight < 0:
        raise ValueError("invalid average-weight parameters")
    if average_weight > m:
        return Fraction(0)

    # At fixed m all U_bal(m,d,W) values share the same prefix sum.  Building
    # it once changes the truncation minimization from quadratic to linear in
    # the number of relevant weight levels.
    prefix = [Fraction(0)] * ((m + d - 1) // d + 1)
    running = Fraction(0)
    for level in range(1, len(prefix)):
        prefix[level] = running
        running += Fraction(d * comb(m, level), comb(d * level, level))

    def capacity_at(cutoff: int) -> Fraction:
        if cutoff <= 0:
            return Fraction(0)
        level = (cutoff + d - 1) // d
        remainder = cutoff - d * (level - 1)
        return prefix[level] + Fraction(
            remainder * comb(m, level), comb(cutoff, level)
        )

    capacity = capacity_at(m)
    for cutoff in range(floor(average_weight), m):
        denominator = cutoff + 1 - average_weight
        if denominator <= 0:
            continue
        candidate = (
            Fraction(cutoff + 1, 1) / denominator
            * capacity_at(cutoff)
        )
        capacity = min(capacity, candidate)
    return capacity


@lru_cache(maxsize=None)
def average_weight_lower_bound(
    n: int, d: int, average_weight: Fraction | int
) -> int:
    """Invert the article's average-column-weight capacity bound."""
    average_weight = Fraction(average_weight)
    if n < 1 or d < 1 or average_weight < 0:
        raise ValueError("invalid average-weight parameters")
    low = max(0, (average_weight.numerator + average_weight.denominator - 1)
              // average_weight.denominator)
    high = max(1, low)
    while average_weight_capacity(high, d, average_weight) < n:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if average_weight_capacity(middle, d, average_weight) >= n:
            high = middle
        else:
            low = middle + 1
    return low


def _binomial_ratio(numerator_n: int, denominator_n: int, k: int) -> float:
    """Return binom(numerator_n, k) / binom(denominator_n, k) stably."""
    if k < 0 or numerator_n < k or denominator_n < numerator_n:
        return 0.0
    if numerator_n == denominator_n:
        return 1.0
    logarithm = (
        lgamma(numerator_n + 1) - lgamma(numerator_n - k + 1)
        - lgamma(denominator_n + 1) + lgamma(denominator_n - k + 1)
    )
    return exp(logarithm)


@lru_cache(maxsize=None)
def average_weight_profile_relaxation_capacity(
    m: int, d: int, average_weight: Fraction | int
) -> float:
    """Upper-bound the size using the article's direct average-profile LP.

    The article requires integral local own-subset profiles.  This routine
    relaxes those local profiles to real vectors.  The resulting LP is still
    a valid code-size upper bound (and hence gives a valid row lower bound),
    while avoiding an infeasible enumeration of all integral profiles.

    Only the two endpoint profiles of the balanced local inequality are
    needed.  Subset capacities are divided by ``binom(m, j)`` for numerical
    stability, so every right-hand side equals one.
    """
    average_weight = Fraction(average_weight)
    if m < 0 or d < 1 or average_weight < 0:
        raise ValueError("invalid average-weight parameters")
    if average_weight > m:
        return 0.0
    if m == 0:
        return 0.0

    # Each column is (weight v, used level j, normalized own-set demand).
    profiles: list[tuple[int, int, float]] = []
    for v in range(1, m + 1):
        t, remainder = divmod(v, d)
        if t > 0 and d - remainder > 0:
            profiles.append((
                v, t,
                _binomial_ratio(v, m, t) / (d - remainder),
            ))
        if remainder > 0:
            level = t + 1
            profiles.append((
                v, level,
                _binomial_ratio(v, m, level) / remainder,
            ))

    # A normalized binomial ratio can underflow only when its reciprocal is
    # beyond the floating-point range.  If such a zero-cost profile does not
    # consume the average-weight budget either, the representable capacity is
    # effectively infinite; passing it to HiGHS would look spuriously
    # unbounded.
    average = float(average_weight)
    if any(demand == 0.0 and v <= average for v, _level, demand in profiles):
        return float("inf")

    row_indices: list[int] = []
    column_indices: list[int] = []
    coefficients: list[float] = []
    for column, (v, level, demand) in enumerate(profiles):
        row_indices.extend((level - 1, m))
        column_indices.extend((column, column))
        coefficients.extend((demand, v - average))
    matrix = coo_matrix(
        (coefficients, (row_indices, column_indices)),
        shape=(m + 1, len(profiles)),
    ).tocsr()
    rhs = np.concatenate((np.ones(m), np.zeros(1)))
    solution = linprog(
        -np.ones(len(profiles)), A_ub=matrix, b_ub=rhs,
        bounds=(0, None), method="highs",
    )
    if solution.status == 3:
        # The mathematical LP is bounded because every profile consumes a
        # positive subset capacity.  HiGHS can nevertheless classify it as
        # unbounded after coefficients below its feasibility tolerance have
        # effectively vanished.  Returning infinity is conservative.
        return float("inf")
    if not solution.success:
        raise RuntimeError(f"average-profile LP failed: {solution.message}")
    return max(0.0, -float(solution.fun))


@lru_cache(maxsize=None)
def average_weight_profile_relaxation_lower_bound(
    n: int, d: int, average_weight: Fraction | int
) -> int:
    """Invert :func:`average_weight_profile_relaxation_capacity` in ``m``."""
    average_weight = Fraction(average_weight)
    if n < 1 or d < 1 or average_weight < 0:
        raise ValueError("invalid average-weight parameters")
    low = max(0, (average_weight.numerator + average_weight.denominator - 1)
              // average_weight.denominator)
    high = max(1, low)
    while average_weight_profile_relaxation_capacity(high, d, average_weight) < n:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if average_weight_profile_relaxation_capacity(middle, d, average_weight) >= n:
            high = middle
        else:
            low = middle + 1
    return low


def combined_average_weight_lower_bound(
    n: int, d: int, average_weight: Fraction | int, general_lower_bound: int
) -> int:
    """Return ``max(general_lower_bound, M_avg)`` without needless inversion."""
    if general_lower_bound < 0:
        raise ValueError("general_lower_bound must be nonnegative")
    average_weight = Fraction(average_weight)
    low = max(
        general_lower_bound,
        (average_weight.numerator + average_weight.denominator - 1)
        // average_weight.denominator,
    )
    def capacity(m: int) -> float:
        return min(
            _average_weight_capacity_heuristic(m, d, average_weight),
            average_weight_profile_relaxation_capacity(m, d, average_weight),
        )

    if capacity(low) >= n:
        return low
    high = max(1, low)
    while capacity(high) < n:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if capacity(middle) >= n:
            high = middle
        else:
            low = middle + 1
    return low


@lru_cache(maxsize=None)
def _average_weight_capacity_heuristic(
    m: int, d: int, average_weight: Fraction | int
) -> float:
    """First-local-minimum heuristic used for the article measurements.

    Starting at the smallest admissible integer cutoff, increase ``W`` while
    the Markov-truncated capacity strictly decreases.  Stop at the first
    non-improving cutoff.  The exhaustive exact implementation remains
    available as :func:`average_weight_capacity` for verification.
    """
    average_weight = Fraction(average_weight)
    if average_weight > m:
        return 0.0
    prefix = [0.0]
    running = 0.0

    def combination_ratio(n1: int, k1: int, n2: int, k2: int) -> float:
        logarithm = (
            lgamma(n1 + 1) - lgamma(k1 + 1) - lgamma(n1 - k1 + 1)
            - lgamma(n2 + 1) + lgamma(k2 + 1) + lgamma(n2 - k2 + 1)
        )
        return float("inf") if logarithm > 709 else exp(logarithm)

    def ensure_prefix(level: int) -> None:
        nonlocal running
        while len(prefix) <= level:
            current = len(prefix)
            prefix.append(running)
            running += d * combination_ratio(m, current, d * current, current)

    def capacity_at(cutoff: int) -> float:
        if cutoff <= 0:
            return 0.0
        level = (cutoff + d - 1) // d
        ensure_prefix(level)
        remainder = cutoff - d * (level - 1)
        tail = remainder * combination_ratio(m, level, cutoff, level)
        return prefix[level] + tail

    average = float(average_weight)
    best = float("inf")
    for cutoff in range(floor(average_weight), m):
        denominator = cutoff + 1 - average
        if denominator > 0:
            candidate = (cutoff + 1) / denominator * capacity_at(cutoff)
            if candidate >= best:
                return best
            best = candidate
    # If no local increase occurred before all possible weights were covered,
    # compare with the direct maximum-weight certificate at W=m.
    return min(best, capacity_at(m))


@lru_cache(maxsize=None)
def recursive_maximum_weight_lower_bound(n: int, d: int, maximum_weight: int) -> int:
    """Combine a supplied weight budget with the balanced DR recursion.

    At a candidate length ``m``, deletion of a maximum-weight column gives
    the additional budget ``m - M_BDR(n-1,d-1)``.  Monotonicity of
    ``balanced_partition_capacity`` allows both restrictions to be imposed
    by taking their minimum.
    """
    if n < 1 or d < 1 or maximum_weight < 1:
        raise ValueError("n, d, and maximum_weight must be positive")
    if d == 1:
        return max(
            sperner_lower_bound(n),
            maximum_weight_lower_bound(n, d, maximum_weight),
        )
    recursive_rows = balanced_dr_lower_bound(n - 1, d - 1)
    m = recursive_rows + 1
    while True:
        effective_weight = min(maximum_weight, m - recursive_rows)
        if balanced_partition_capacity(m, d, effective_weight) >= n:
            return m
        m += 1


@lru_cache(maxsize=None)
def sperner_lower_bound(n: int) -> int:
    """Minimum row count for an antichain of at least ``n`` columns."""
    if n < 1:
        raise ValueError("n must be positive")
    m = 0
    while comb(m, m // 2) < n:
        m += 1
    return m


@lru_cache(maxsize=None)
def balanced_dr_lower_bound(n: int, d: int) -> int:
    """Compute the recursive balanced D'yachkov-Rykov bound M_BDR(n,d)."""
    if n < 1 or d < 1:
        raise ValueError("n and d must be positive")
    if n == 1:
        return 1
    if d == 1:
        return sperner_lower_bound(n)

    recursive_rows = balanced_dr_lower_bound(n - 1, d - 1)
    # m must leave at least recursive_rows coordinates after deleting the
    # maximum-weight column. Start at the same unavoidable row count.
    m = recursive_rows + 1
    while True:
        maximum_weight = m - recursive_rows
        if balanced_partition_capacity(m, d, maximum_weight) >= n:
            return m
        m += 1


def balanced_dr_certificate(n: int, d: int) -> dict[str, int | str]:
    """Return the bound and the parameters witnessing its final comparison."""
    m = balanced_dr_lower_bound(n, d)
    if d == 1:
        return {"n": n, "d": d, "lower_bound": m, "base": "Sperner"}
    recursive_rows = balanced_dr_lower_bound(n - 1, d - 1)
    maximum_weight = m - recursive_rows
    capacity = balanced_partition_capacity(m, d, maximum_weight)
    return {
        "n": n,
        "d": d,
        "lower_bound": m,
        "recursive_lower_bound": recursive_rows,
        "maximum_weight": maximum_weight,
        "capacity_numerator": capacity.numerator,
        "capacity_denominator": capacity.denominator,
    }
