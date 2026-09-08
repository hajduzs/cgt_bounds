from fractions import Fraction

import pytest

from cgt.bounds.balanced_dyachkov_rykov import (
    balanced_dr_certificate,
    balanced_dr_lower_bound,
    balanced_partition_capacity,
    maximum_weight_lower_bound,
    sperner_lower_bound,
)


def test_balanced_capacity_uses_exact_rationals():
    assert balanced_partition_capacity(7, 3, 4) == Fraction(21, 2)


@pytest.mark.parametrize("n, expected", [(1, 0), (2, 2), (3, 3), (6, 4), (10, 5)])
def test_sperner_base(n, expected):
    assert sperner_lower_bound(n) == expected


def test_maximum_weight_inversion():
    assert maximum_weight_lower_bound(100, 2, 5) == 16
    assert maximum_weight_lower_bound(1000, 2, 7) == 28


def test_balanced_dr_is_recursive_and_certified():
    value = balanced_dr_lower_bound(100, 2)
    certificate = balanced_dr_certificate(100, 2)
    assert value == certificate["lower_bound"]
    assert certificate["recursive_lower_bound"] == sperner_lower_bound(99)
    assert balanced_partition_capacity(value, 2, certificate["maximum_weight"]) >= 100


def test_invalid_parameters():
    with pytest.raises(ValueError):
        balanced_dr_lower_bound(0, 2)
