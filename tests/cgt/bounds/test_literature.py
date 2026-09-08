from fractions import Fraction

import pytest

from cgt.bounds.literature import (
    inan_kairouz_ozgur_capacity,
    inan_kairouz_ozgur_lower_bound,
    shangguan_ge_constant_weight_lower_bound,
    shangguan_ge_lower_bound,
    shangguan_ge_maximum_weight_lower_bound,
    shangguan_ge_weight_aware_lower_bound,
)


def test_shangguan_ge_uses_exact_irrational_ceiling():
    # ceil(((15 + sqrt(33))/24) * d^2)
    assert shangguan_ge_lower_bound(10_000, 2) == 4
    assert shangguan_ge_lower_bound(10_000, 10) == 87
    assert shangguan_ge_lower_bound(50, 10) == 50


def test_shangguan_ge_weight_dependent_corollaries():
    assert shangguan_ge_constant_weight_lower_bound(1000, 3, 4) == 16
    assert shangguan_ge_constant_weight_lower_bound(1000, 3, 5) == 0
    assert shangguan_ge_maximum_weight_lower_bound(1000, 3, 5) == 14
    assert shangguan_ge_maximum_weight_lower_bound(1000, 3, 6) == 0
    assert shangguan_ge_weight_aware_lower_bound(1000, 3, 4) == 16
    assert shangguan_ge_weight_aware_lower_bound(10, 3, 4) == 10


def test_shangguan_ge_weight_aware_theorem_1_1_case():
    # Theorem 1.1: constant column weight d+1 forces t >= (d+1)^2.
    assert shangguan_ge_weight_aware_lower_bound(10_000, 10, 11) == 121
    assert shangguan_ge_weight_aware_lower_bound(50, 10, 11) == 50


def test_shangguan_ge_weight_aware_corollary_3_5_case():
    # Corollary 3.5 (integer form): weight <= floor(5d/3) forces
    # t >= d^2+d+2, i.e. t > d^2+d+1.
    assert shangguan_ge_weight_aware_lower_bound(10_000, 10, 16) == 112
    assert shangguan_ge_weight_aware_lower_bound(50, 10, 16) == 50
    # d=10: floor(5*10/3) = 16, so weight 17 no longer qualifies.
    assert shangguan_ge_weight_aware_lower_bound(10_000, 10, 17) == 0


def test_shangguan_ge_weight_aware_takes_the_larger_applicable_term():
    # d=2: d+1=3 triggers Theorem 1.1 ((d+1)^2=9); floor(5*2/3)=3 also
    # triggers Corollary 3.5 (d^2+d+2=8). Theorem 1.1's term wins.
    assert shangguan_ge_weight_aware_lower_bound(10_000, 2, 3) == 9


@pytest.mark.parametrize("args", [(0, 2, 3), (10, 0, 3), (10, 2, 0)])
def test_shangguan_ge_weight_aware_rejects_invalid_parameters(args):
    with pytest.raises(ValueError):
        shangguan_ge_weight_aware_lower_bound(*args)


def test_iko_identity_and_private_pair_cases():
    assert inan_kairouz_ozgur_lower_bound(100, 2, 2) == 100
    assert inan_kairouz_ozgur_lower_bound(100, 2, 3) == 25
    assert inan_kairouz_ozgur_lower_bound(100, 3, 4) == 35


def test_iko_general_capacity_is_exact_fraction():
    # ell=2: C(m,1)+C(m,2)+C(m,3)/C(d+2,3)
    expected = Fraction(10 + 45) + Fraction(120, 4)
    assert inan_kairouz_ozgur_capacity(10, 2, 5) == expected
    lower = inan_kairouz_ozgur_lower_bound(86, 2, 5)
    assert lower == 11
    assert inan_kairouz_ozgur_capacity(lower - 1, 2, 5) < 86
    assert inan_kairouz_ozgur_capacity(lower, 2, 5) >= 86


@pytest.mark.parametrize("args", [(0, 2, 3), (10, 0, 3), (10, 2, 0)])
def test_iko_rejects_invalid_parameters(args):
    with pytest.raises(ValueError):
        inan_kairouz_ozgur_lower_bound(*args)
