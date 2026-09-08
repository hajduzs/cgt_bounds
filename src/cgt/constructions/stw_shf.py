"""Stinson--Trung--Wei recursive separating-hash construction.

Reference: D. R. Stinson, T. van Trung & R. Wei (2000). "Secure
Frameproof Codes, Key Distribution Patterns, Group Testing Algorithms
and Related Structures." Journal of Statistical Planning and
Inference, 86(2):595-617. DOI: 10.1016/S0378-3758(99)00131-7.

For type ``{1, d}``, Theorem 4.8 of Stinson, Trung and Wei turns an
``SHF(N0; n0, q, {1,d})`` into

    SHF((d + 1)^j N0; n0^(2^j), q, {1,d})

provided ``gcd(n0, d!) = 1``.  One-hot expansion of the alphabet gives a
binary d-disjunct matrix with ``m = q N`` tests and column weight ``N``.
"""

from __future__ import annotations

from math import factorial, gcd

from cgt.constructions.ks64 import SuperimposedCode, _qary_mds_words
from cgt.tables.ks import _is_prime_power


def _recursive_step(
    words: list[tuple[int, ...]], d: int
) -> list[tuple[int, ...]]:
    """Apply the column-permutation array in STW Theorem 4.8 once."""
    n0 = len(words)
    result: list[tuple[int, ...]] = []
    # Columns are indexed by (block, local).  In row block i, block j is
    # the base array after the cyclic column shift i*j (mod n0).
    for block in range(n0):
        for local in range(n0):
            result.append(
                tuple(
                    symbol
                    for slope in range(d + 1)
                    for symbol in words[(local - slope * block) % n0]
                )
            )
    return result


def stw_recursive_shf(
    q: int,
    dimension: int,
    d: int,
    levels: int,
    n: int | None = None,
) -> SuperimposedCode:
    """Materialize the RS-seeded STW construction and one-hot expand it.

    ``n`` may truncate the resulting family; deletion of columns preserves
    disjunctness and avoids materializing an unnecessarily wide binary matrix.
    The q-ary recursion itself still has quadratic growth at every level, so
    explicit generation is intended for moderate instances.
    """
    if q < 2 or not _is_prime_power(q):
        raise ValueError("q must be a prime power")
    if dimension < 1 or d < 1 or levels < 0:
        raise ValueError("dimension and d must be positive and levels nonnegative")
    base_length = 1 + d * (dimension - 1)
    if base_length > q:
        raise ValueError("the RS seed requires 1 + d*(dimension-1) <= q")
    base_size = q**dimension
    if gcd(base_size, factorial(d)) != 1:
        raise ValueError("STW requires gcd(q**dimension, d!) = 1")

    words = _qary_mds_words(q, dimension, base_length)
    for _ in range(levels):
        words = _recursive_step(words, d)
    if n is not None:
        if n < 1 or n > len(words):
            raise ValueError(f"n must lie between 1 and {len(words)}")
        words = words[:n]

    columns = [
        frozenset(row * q + symbol for row, symbol in enumerate(word))
        for word in words
    ]
    weight = base_length * (d + 1) ** levels
    return SuperimposedCode(
        columns,
        q * weight,
        d,
        "disjunct",
        "stinson-trung-wei-recursive-shf",
        {
            "q": q,
            "dimension": dimension,
            "levels": levels,
            "base_N": base_length,
            "base_n": base_size,
            "full_n": base_size ** (2**levels),
            "theorem": "Stinson--Trung--Wei 2000, Theorem 4.8",
        },
    )


def stw_parameters(
    n: int, d: int, maximum_q: int = 257, minimum_levels: int = 1
) -> tuple[int, int, int, int, int, int]:
    """Return the best ``(m,w,q,k,j,capacity)`` RS-seeded STW candidate.

    The default requires at least one genuinely recursive step.  Setting
    ``minimum_levels=0`` also admits the unmodified RS seed.
    """
    if n < 1 or d < 1 or maximum_q < 2 or minimum_levels < 0:
        raise ValueError("n, d, maximum_q, and minimum_levels must be nonnegative")
    candidates = []
    for q in range(2, maximum_q + 1):
        if not _is_prime_power(q) or gcd(q, factorial(d)) != 1:
            continue
        dimension = 1
        while 1 + d * (dimension - 1) <= q:
            base_N = 1 + d * (dimension - 1)
            capacity = q**dimension
            levels = 0
            weight = base_N
            while capacity < n or levels < minimum_levels:
                capacity *= capacity
                weight *= d + 1
                levels += 1
            candidates.append((q * weight, weight, q, dimension, levels, capacity))
            dimension += 1
    if not candidates:
        raise RuntimeError(f"no STW parameter found up to q={maximum_q}")
    return min(candidates)
