import pytest

from cgt.constructions.ks64 import is_d_disjunct
from cgt.constructions.stw_shf import stw_parameters, stw_recursive_shf


def test_stw_one_level_is_explicitly_two_disjunct():
    code = stw_recursive_shf(q=3, dimension=1, d=2, levels=1)
    assert (code.m, code.N, code.w) == (9, 9, 3)
    assert is_d_disjunct(code)


def test_stw_can_truncate_columns():
    code = stw_recursive_shf(q=3, dimension=1, d=2, levels=1, n=7)
    assert code.N == 7
    assert is_d_disjunct(code)


def test_stw_enforces_coprimality_condition():
    with pytest.raises(ValueError, match="gcd"):
        stw_recursive_shf(q=4, dimension=1, d=2, levels=1)


def test_stw_parameter_search_finds_recursive_improvement():
    assert stw_parameters(10_000, 2)[:5] == (75, 15, 5, 3, 1)
