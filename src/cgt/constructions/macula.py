"""Macula's inclusion-matrix construction.

Reference: A. J. Macula (1996). "A Simple Construction of d-Disjunct
Matrices with Certain Constant Weights." Discrete Mathematics,
162(1-3):311-312. DOI: 10.1016/0012-365X(95)00296-9.
"""

from __future__ import annotations

from math import comb


def macula_parameters(n: int, d: int) -> tuple[int, int, int, int]:
    """Return the shortest Macula inclusion construction ``(m,w,v,k)``.

    Rows are the d-subsets of a v-set and columns are its k-subsets.  Thus
    m=C(v,d), N=C(v,k), w=C(k,d), with d<k<v.
    """
    if n < 2 or d < 1:
        raise ValueError("require n >= 2 and d >= 1")
    v = d + 2
    while True:
        choices = range(d + 1, v)
        k = max(choices, key=lambda value: comb(v, value))
        if comb(v, k) >= n:
            return comb(v, d), comb(k, d), v, k
        v += 1
