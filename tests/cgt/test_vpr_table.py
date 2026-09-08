import pytest

from cgt.tables.vpr import vpr_parameters, vpr_rows


def test_vpr_reproduces_survey_column():
    rows = vpr_rows()
    assert [row.m for row in rows] == [105, 140, 182, 231, 273, 364, 455, 910, 1358]
    assert [row.w for row in rows] == [21, 20, 26, 33, 39, 52, 65, 130, 194]


def test_vpr_balanced_partition_gaps():
    rows = vpr_rows()
    assert [row.lower_bound for row in rows] == [13, 22, 32, 42, 52, 71, 92, 193, 295]
    assert rows[-1].gap_ratio == pytest.approx(1358 / 295)


def test_vpr_d3_reproduces_survey_column():
    rows = vpr_rows(d=3)
    assert [row.m for row in rows] == [287, 451, 682, 902, 1122, 2244, 3366]
    assert [row.q for row in rows] == [7, 11, 11, 11, 11, 11, 11]
    assert [row.lower_bound for row in rows] == [18, 48, 79, 111, 143, 307, 473]


def test_vpr_rejects_unsupported_d():
    with pytest.raises(ValueError, match="d=2 and d=3"):
        vpr_rows(d=5)


def test_vpr_single_parameter_choice():
    assert vpr_parameters(10**4, 2, 7) == (182, 26)


def test_vpr_rejects_relative_distance_outside_gv_range():
    with pytest.raises(ValueError, match="GV range"):
        vpr_parameters(100, 10, 2)
