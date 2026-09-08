import pytest

from cgt.tables.ks import ks_b_parameters, ks_c_parameters, ks_rows


@pytest.mark.parametrize(
    "n, expected",
    [
        (100, (25, 5)), (10**3, (49, 7)), (10**4, (77, 11)),
        (10**5, (99, 11)), (10**6, (121, 11)), (10**8, (208, 16)),
        (10**10, (285, 19)), (10**20, (729, 27)), (10**30, (1517, 41)),
    ],
)
def test_ks_b_reproduces_survey(n, expected):
    assert ks_b_parameters(n, 2) == expected


def test_ks_rows_reproduce_reference_columns():
    rows = ks_rows()
    assert [row.ks_c_rs_m for row in rows] == [21, 49, 75, 99, 108, 156, 225, 475, 725]
    assert [row.ks_c_m for row in rows] == [21, 39, 63, 81, 99, 143, 180, 368, 555]
    assert [row.ks_c_weight for row in rows] == [5, 15, 21, 27, 33, 39, 45, 115, 185]
    assert [row.as_dict()["ks_c_lb"] for row in rows] == [16, 23, 32, 42, 52, 72, 92, 193, 295]


def test_unknown_reference_row_is_explicit():
    with pytest.raises(ValueError, match="no survey KS-C reference"):
        ks_rows(ns=[101])


def test_ks_c_fills_the_intervening_d2_decades():
    assert ks_c_parameters(10**7, 2) == (117, 39)
    assert ks_c_parameters(10**9, 2) == (165, 45)


def test_ks_c_searches_stored_inner_codes_beyond_d2():
    # The search retains the one-hot specialization as a certified fallback.
    assert ks_c_parameters(10**3, 3)[0] <= ks_b_parameters(10**3, 3)[0]
