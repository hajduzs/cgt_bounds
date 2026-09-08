from fractions import Fraction
import pytest

from cgt.bounds.classical import (
    c_dr_capacity,
    dyachkov_rykov_capacity,
    dyachkov_rykov_lower_bound,
    kautz_singleton_lower_bound,
    maximum_weight_dyachkov_rykov_lower_bound,
)


def test_exact_u_dr_capacity_and_closed_relaxation_are_distinct():
    # U_DR(10,2,5) = 10 + C(10,2)/C(2,1)
    #                        + C(10,2)/C(3,1)
    #                        + C(10,3)/C(4,2).
    assert dyachkov_rykov_capacity(10, 2, 5) == Fraction(135, 2)
    assert dyachkov_rykov_capacity(10, 2, 5) <= c_dr_capacity(10, 2, 5)


def test_kautz_singleton_reproduces_the_counting_values_for_d2():
    assert [kautz_singleton_lower_bound(10**e, 2) for e in range(2, 11)] == [
        13, 19, 26, 33, 39, 46, 53, 59, 66,
    ]


def test_classical_dr_is_no_stronger_than_balanced_dr():
    from cgt.bounds.balanced_dyachkov_rykov import balanced_dr_lower_bound

    for d in (2, 3, 5, 10):
        for exponent in (2, 5):
            n = 10**exponent
            assert dyachkov_rykov_lower_bound(n, d) <= balanced_dr_lower_bound(n, d)


def test_maximum_weight_dr_bound_inverts_its_capacity():
    lower = maximum_weight_dyachkov_rykov_lower_bound(1000, 2, 7)
    assert dyachkov_rykov_capacity(lower, 2, 7) >= 1000
    assert dyachkov_rykov_capacity(lower - 1, 2, 7) < 1000


def test_recursive_weight_bound_dominates_its_two_components():
    # Cost grows ~10x per exponent step (n=10**10 alone takes ~4s per call,
    # vs ~0.0004s at n=10**2), so keep exponents small here -- (2, 5) is
    # cheap and still exercises both d extremes. See the _full_sweep
    # companion below (slow) for the exhaustive check up to n=10**10.
    from cgt.bounds.balanced_dyachkov_rykov import (
        balanced_dr_lower_bound,
        maximum_weight_lower_bound,
        recursive_maximum_weight_lower_bound,
    )

    for d in (2, 10):
        for exponent in (2, 5):
            n = 10**exponent
            for weight in (d + 1, 2 * d + 3):
                combined = recursive_maximum_weight_lower_bound(n, d, weight)
                assert combined >= balanced_dr_lower_bound(n, d)
                assert combined >= maximum_weight_lower_bound(n, d, weight)


@pytest.mark.slow
def test_recursive_weight_bound_dominates_its_two_components_full_sweep():
    from cgt.bounds.balanced_dyachkov_rykov import (
        balanced_dr_lower_bound,
        maximum_weight_lower_bound,
        recursive_maximum_weight_lower_bound,
    )

    for d in (2, 3, 5, 10):
        for exponent in range(2, 11):
            n = 10**exponent
            for weight in (d + 1, 2 * d + 3):
                combined = recursive_maximum_weight_lower_bound(n, d, weight)
                assert combined >= balanced_dr_lower_bound(n, d)
                assert combined >= maximum_weight_lower_bound(n, d, weight)
