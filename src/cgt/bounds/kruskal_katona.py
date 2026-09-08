"""Finite Kruskal-Katona refinement.

Reference:
    - J. B. Kruskal (1963). "The Number of Simplices in a Complex." In
      Mathematical Optimization Techniques (R. Bellman, ed.), University
      of California Press, 251-278. DOI: 10.1525/9780520319875-014.
      (the shadow-minimization theorem this module's LP is built on)
    - Refinement to a finite disjunct-matrix bound: this paper (see the
      accompanying article, Appendix D).
"""

from __future__ import annotations

from functools import lru_cache
from math import ceil, comb

import numpy as np
from scipy.optimize import linprog

from .balanced_dyachkov_rykov import (
    balanced_partition_capacity,
    maximum_weight_lower_bound,
)


PROFILE_LIMIT = 5_000


def kk_upper_shadow(v: int, j: int, p: int) -> int:
    """Exact minimum upper shadow of ``p`` j-subsets of a v-set."""
    if p <= 0:
        return 0
    if not 0 <= j < v or p > comb(v, j):
        raise ValueError("invalid Kruskal--Katona parameters")
    k = v - j
    remainder = p
    previous = v + 1
    shadow = 0
    for rank in range(k, 0, -1):
        a = min(previous - 1, v)
        while a >= rank and comb(a, rank) > remainder:
            a -= 1
        if a >= rank:
            remainder -= comb(a, rank)
            shadow += comb(a, rank - 1)
            previous = a
        if remainder == 0:
            break
    return shadow


def _minimal_profiles(d: int, v: int) -> tuple[tuple[int, ...], ...] | None:
    """Enumerate nondominated profiles needed by the article's finite LP."""
    t, r = divmod(v, d)
    if t == 0:
        profile = [0] * (v + 1)
        profile[1] = 1
        for j in range(1, v):
            profile[j + 1] = kk_upper_shadow(v, j, profile[j])
        return (tuple(profile),)

    count = comb(v, t) + 1
    if count > PROFILE_LIMIT:
        return None
    ct = comb(v, t)
    ct1 = comb(v, t + 1)
    profiles = []
    for qt in range(ct + 1):
        required = 0
        if r:
            required = max(0, ceil((1 - (d - r) * qt / ct) * ct1 / r - 1e-12))
        elif qt < ct / d:
            continue
        qt1 = max(kk_upper_shadow(v, t, qt), required)
        if qt1 > ct1:
            continue
        q = [0] * (v + 1)
        q[t], q[t + 1] = qt, qt1
        for j in range(t + 1, v):
            q[j + 1] = kk_upper_shadow(v, j, q[j])
        if q[v] == 1:
            profiles.append(tuple(q))
    # For fixed q_t, choosing any larger later coordinate is dominated by
    # the minimally propagated profile above.  Retaining this one-dimensional
    # frontier avoids an otherwise quadratic dominance pass.
    return tuple(dict.fromkeys(profiles))


@lru_cache(maxsize=None)
def kk_profiles(d: int, maximum_weight: int):
    profiles = []
    for v in range(1, maximum_weight + 1):
        local = _minimal_profiles(d, v)
        if local is None or len(profiles) + len(local) > PROFILE_LIMIT:
            return None
        profiles.extend(local)
    return tuple(profiles)


@lru_cache(maxsize=None)
def kk_exact_weight_limit(d: int) -> int:
    """Largest weight for which the configured exact profile LP is feasible."""
    profiles = 0
    weight = 1
    while True:
        local = _minimal_profiles(d, weight)
        if local is None or profiles + len(local) > PROFILE_LIMIT:
            return weight - 1
        profiles += len(local)
        weight += 1


@lru_cache(maxsize=None)
def kk_capacity(m: int, d: int, maximum_weight: int) -> float:
    """Return U_KK, falling back to U_bal when exact enumeration is too large."""
    maximum_weight = min(m, maximum_weight)
    profiles = kk_profiles(d, maximum_weight)
    if profiles is None:
        return float(balanced_partition_capacity(m, d, maximum_weight))
    matrix = np.array(
        [[profile[j] if j < len(profile) else 0 for profile in profiles]
         for j in range(1, maximum_weight + 1)],
        dtype=float,
    )
    result = linprog(
        -np.ones(len(profiles)), A_ub=matrix,
        b_ub=np.array([comb(m, j) for j in range(1, maximum_weight + 1)], dtype=float),
        # q_v=1 supplies this bound implicitly.  Passing it explicitly avoids
        # false "unbounded" reports from HiGHS when binomial capacities span
        # many orders of magnitude.
        bounds=[(0, float(comb(m, len(profile) - 1))) for profile in profiles],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(f"Kruskal--Katona LP failed: {result.message}")
    return float(-result.fun)


def kk_maximum_weight_lower_bound(n: int, d: int, maximum_weight: int) -> int:
    """Invert the article's multilevel Kruskal--Katona certificate."""
    if n < 1 or d < 1 or maximum_weight < 1:
        raise ValueError("n, d, and maximum_weight must be positive")
    # U_KK <= U_bal, hence its inverse cannot start below the balanced bound.
    m = maximum_weight_lower_bound(n, d, maximum_weight)
    while kk_capacity(m, d, maximum_weight) + 1e-8 < n:
        m += 1
    return m


@lru_cache(maxsize=None)
def kk_unrestricted_lower_bound(n: int, d: int) -> int:
    """Invert the KK certificate without a supplied weight restriction.

    At a candidate length ``m`` every binary column has weight at most
    ``m``.  Thus ``kk_capacity(m,d,m)`` is the appropriate unrestricted
    capacity.  This is the KK component of the article's ``LB_infty``.
    """
    if n < 1 or d < 1:
        raise ValueError("n and d must be positive")
    exact_limit = kk_exact_weight_limit(d)

    def capacity(m: int) -> float:
        if m <= exact_limit:
            return kk_capacity(m, d, m)
        return float(balanced_partition_capacity(m, d, m))

    m = 1
    while capacity(m) + 1e-8 < n:
        m += 1
    return m


def kk_was_exact(d: int, maximum_weight: int) -> bool:
    """Whether the profile LP, rather than U_bal fallback, was evaluated."""
    return kk_profiles(d, maximum_weight) is not None
