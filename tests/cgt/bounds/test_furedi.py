import pytest

from cgt.bounds.furedi import furedi_capacity, furedi_certificate, furedi_lower_bound


def test_d2_capacity_formula():
    assert furedi_capacity(21, 2) == 2 + 116280


def test_strict_inequality_is_inverted_exactly():
    bound = furedi_lower_bound(100000, 2)
    assert bound == 21
    assert 100000 >= furedi_capacity(bound - 1, 2)
    assert 100000 < furedi_capacity(bound, 2)


def test_certificate():
    assert furedi_certificate(1000, 2) == {
        "name": "Furedi 1996 r-cover-free bound",
        "n": 1000,
        "d": 2,
        "lower_bound": 14,
        "binomial_layer": 4,
        "strict_capacity": 1003,
    }


def test_invalid_parameters():
    with pytest.raises(ValueError):
        furedi_lower_bound(0, 2)
