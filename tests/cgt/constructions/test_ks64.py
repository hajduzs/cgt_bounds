import networkx as nx
import pytest

from cgt.constructions.ks64 import (
    bch_parity_ud2,
    bch_parity_ud3,
    graph_ud2,
    hamming_weight3_zfd,
    identity_zfd,
    is_d_disjunct,
    is_d_separable,
    iterated_qary_composition,
    ks64_ud2_seed,
    latin_square_zfd,
    moore_graph_ud2,
    projective_plane_zfd,
    qary_mds_zfd,
    split_field_bibd_ud2,
    split_projective_plane_ud2,
)


@pytest.mark.parametrize(
    "code, expected",
    [
        (identity_zfd(5, 2), (5, 5, 1)),
        (hamming_weight3_zfd(3), (7, 7, 3)),
        (qary_mds_zfd(3, 2, 2), (9, 9, 3)),
        (latin_square_zfd(4), (12, 16, 3)),
        (iterated_qary_composition(3, 2, 2, 2), (27, 81, 9)),
        (projective_plane_zfd(2), (7, 7, 3)),
    ],
)
def test_explicit_disjunct_constructions(code, expected):
    assert (code.m, code.N, code.w) == expected
    assert code.guarantee == "disjunct"
    assert code.historical_name == f"ZFD_{code.d}"
    assert is_d_disjunct(code)


@pytest.mark.parametrize(
    "code, expected",
    [
        (bch_parity_ud2(3), (12, 7)),
        (bch_parity_ud3(3), (144, 7)),
        (moore_graph_ud2(2), (5, 5)),
        (moore_graph_ud2(3), (10, 15)),
        (ks64_ud2_seed(), (7, 9)),
        (split_projective_plane_ud2(2), (14, 21)),
    ],
)
def test_explicit_separable_constructions(code, expected):
    assert (code.m, code.N) == expected
    assert code.guarantee == "separable"
    assert code.historical_name == f"UD_{code.d}"
    assert is_d_separable(code)


def test_graph_construction_rejects_short_cycles():
    with pytest.raises(ValueError, match="girth"):
        graph_ud2(nx.cycle_graph(4))


def test_split_field_bibd_construction():
    fano_blocks = projective_plane_zfd(2).columns
    code = split_field_bibd_ud2(7, fano_blocks)
    assert (code.m, code.N, code.w) == (14, 21, 2)
    assert is_d_separable(code)
