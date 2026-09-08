from math import comb, e

import pytest

from cgt.constructions.gargano import gargano_length, gargano_low_weight_zfd, gargano_parameters
from cgt.constructions.ks64 import is_d_disjunct


def _theorem_lhs(m: int, n: int, d: int, weight: int) -> float:
    k = d + 1
    ratio = ((weight * d - (weight - 1) / 2) / (m - (weight - 1) / 2)) ** weight
    dependency = k * (comb(n, k) - comb(n - k + 1, k))
    return e * ratio * dependency


@pytest.mark.parametrize("n,d,weight", [(4, 1, 1), (8, 2, 2), (12, 3, 4)])
def test_gargano_length_is_the_least_certified_integer(n, d, weight):
    m = gargano_length(n, d, weight)
    assert _theorem_lhs(m, n, d, weight) <= 1
    assert m == weight or _theorem_lhs(m - 1, n, d, weight) > 1


def test_gargano_resampling_generates_constant_weight_disjunct_code():
    code = gargano_low_weight_zfd(n=6, d=2, weight=2, seed=7)
    assert (code.N, code.w, code.m) == (6, 2, gargano_length(6, 2, 2))
    assert is_d_disjunct(code)


def test_gargano_parameter_search_returns_best_weight_in_requested_range():
    assert gargano_parameters(20, 2, maximum_weight=8) == min(
        (gargano_length(20, 2, weight), weight) for weight in range(1, 9)
    )
