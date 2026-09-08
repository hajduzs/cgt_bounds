from math import comb

from cgt.constructions.macula import macula_parameters


def test_macula_parameters_cover_requested_columns():
    m, w, v, k = macula_parameters(100, 2)
    assert comb(v, k) >= 100
    assert m == comb(v, 2)
    assert w == comb(k, 2)
    assert 2 < k < v
