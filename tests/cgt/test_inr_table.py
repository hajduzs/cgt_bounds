import pytest

from cgt.tables.inr import inr_parameters, inr_rows


def test_inr_single_parameter_choice():
    assert inr_parameters(10**4, 2) == 255125


def test_inr_rows_default_survey():
    rows = inr_rows()
    assert [row.n for row in rows] == [100, 10**3, 10**4, 10**5, 10**6, 10**8, 10**10, 10**20, 10**30]
    assert [row.m for row in rows][:3] == [127563, 191344, 255125]
    assert all(row.m == inr_parameters(row.n, row.d) for row in rows)


def test_inr_rows_gap_columns():
    rows = inr_rows()
    row = rows[0]
    assert row.additive_gap == row.m - row.lower_bound
    assert row.gap_ratio == pytest.approx(row.m / row.lower_bound)


def test_inr_rows_scale_with_d():
    d2 = inr_parameters(10**6, 2)
    d3 = inr_parameters(10**6, 3)
    assert d3 == pytest.approx(d2 * 9 / 4, rel=1e-4)


def test_inr_rejects_invalid_parameters():
    with pytest.raises(ValueError, match="n >= 2 and d >= 2"):
        inr_parameters(1, 2)
    with pytest.raises(ValueError, match="n >= 2 and d >= 2"):
        inr_parameters(100, 1)
