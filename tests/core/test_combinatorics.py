import math

import pytest

from core import binom, binom_float, johnson_bound, johnson_bound_without_floor


def test_binom_matches_math_comb():
    assert binom(10, 3) == math.comb(10, 3)
    assert binom(52, 5) == math.comb(52, 5)


def test_binom_out_of_range_is_zero():
    assert binom(5, -1) == 0
    assert binom(5, 6) == 0
    assert binom(-1, 0) == 0


def test_binom_float_matches_exact_for_small_values():
    assert binom_float(10, 3) == pytest.approx(math.comb(10, 3))


def test_binom_float_out_of_range_is_zero():
    assert binom_float(5, -1) == 0.0
    assert binom_float(5, 6) == 0.0


def test_johnson_bound_known_value():
    # A(20, 4, 5) Johnson upper bound
    assert johnson_bound(20, 4, 5) == 912


def test_johnson_bound_without_floor_final_floor_matches_no_floor_result():
    # johnson_bound_without_floor floors only the final product; johnson_bound
    # floors at every step of the recursion, so the two are not expected to
    # agree with each other -- only within johnson_bound_without_floor itself.
    m, D, w = 20, 4, 5
    exact = johnson_bound_without_floor(m, D, w, no_floor=True)
    assert johnson_bound_without_floor(m, D, w) == math.floor(exact)
