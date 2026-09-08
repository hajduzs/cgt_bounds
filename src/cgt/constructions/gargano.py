"""Low-weight superimposed codes of Gargano, Rescigno, and Vaccaro.

Reference: L. Gargano, A. A. Rescigno & U. Vaccaro (2020). "Low-Weight
Superimposed Codes and Related Combinatorial Structures: Bounds and
Applications." Theoretical Computer Science, 806:655-672.
DOI: 10.1016/j.tcs.2019.10.032.
"""

from __future__ import annotations

from itertools import combinations
from math import ceil, comb, log
from random import Random

from cgt.constructions.ks64 import SuperimposedCode


def gargano_length(n: int, d: int, weight: int) -> int:
    """Return the least integer length certified by Gargano et al., Thm. 1.

    Their ``(k,w,n)`` convention corresponds to ``d=k-1`` in the
    group-testing convention used by this package.
    """
    if n < 2 or d < 1 or weight < 1 or d >= n:
        raise ValueError("require n >= 2, 1 <= d < n, and weight >= 1")
    k = d + 1
    dependency = k * (comb(n, k) - comb(n - k + 1, k))
    if dependency <= 0:
        raise ValueError("the theorem requires a positive dependency count")
    numerator = weight * d - (weight - 1) / 2
    offset = (weight - 1) / 2

    def certified(m: int) -> bool:
        denominator = m - offset
        if m < weight or denominator <= 0:
            return False
        log_lhs = 1 + weight * (log(numerator) - log(denominator)) + log(dependency)
        return log_lhs <= 0

    low = weight - 1
    high = max(weight, 2 * weight)
    while not certified(high):
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if certified(middle):
            high = middle
        else:
            low = middle
    return high


def gargano_parameters(
    n: int, d: int, maximum_weight: int | None = None
) -> tuple[int, int]:
    """Return the best ``(m, weight)`` found from the finite theorem.

    By default the search includes four times the paper's natural
    ``weight=(d+1) log(n)`` scale.  Supplying ``maximum_weight`` makes the
    finite search range explicit and reproducible.
    """
    if n < 2 or d < 1 or d >= n:
        raise ValueError("require n >= 2 and 1 <= d < n")
    limit = (
        ceil(4 * (d + 1) * log(n))
        if maximum_weight is None
        else maximum_weight
    )
    if limit < 1:
        raise ValueError("maximum_weight must be positive")
    return min((gargano_length(n, d, weight), weight) for weight in range(1, limit + 1))


def _first_bad_event(columns: list[frozenset[int]], d: int):
    """Return one violated disjunctness event, or ``None`` if there is none."""
    indices = range(len(columns))
    for target in indices:
        others = [index for index in indices if index != target]
        for covering in combinations(others, d):
            union = frozenset().union(*(columns[index] for index in covering))
            if columns[target] <= union:
                return (target, *covering)
    return None


def gargano_low_weight_zfd(
    n: int,
    d: int,
    weight: int,
    seed: int | None = None,
    max_resamples: int = 100_000,
) -> SuperimposedCode:
    """Generate the constant-weight code using Moser--Tardos resampling.

    Matrix generation is intended for moderate parameters because detecting
    a bad event enumerates ``d``-tuples of columns.
    """
    if max_resamples < 0:
        raise ValueError("max_resamples must be nonnegative")
    m = gargano_length(n, d, weight)
    rng = Random(seed)

    def new_column() -> frozenset[int]:
        return frozenset(rng.sample(range(m), weight))

    columns = [new_column() for _ in range(n)]
    resamples = 0
    while (event := _first_bad_event(columns, d)) is not None:
        if resamples >= max_resamples:
            raise RuntimeError(
                "Gargano LLL resampling did not converge within max_resamples; "
                "increase the limit or choose a different seed"
            )
        for index in set(event):
            columns[index] = new_column()
        resamples += 1

    return SuperimposedCode(
        columns=columns,
        m=m,
        d=d,
        guarantee="disjunct",
        construction="gargano2020-low-weight",
        metadata={
            "weight": weight,
            "seed": seed,
            "resamples": resamples,
            "theorem": "Gargano--Rescigno--Vaccaro 2020, Theorem 1",
        },
    )
