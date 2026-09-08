from cgt.bounds.ruszinko import (
    ruszinko_capacity,
    ruszinko_lower_bound,
    ruszinko_maximum_weight_capacity,
    ruszinko_maximum_weight_lower_bound,
)


def test_ruszinko_bound_inverts_exact_capacity():
    for d in (2, 3, 5, 10):
        lower = ruszinko_lower_bound(1000, d)
        assert ruszinko_capacity(lower, d) >= 1000
        if lower:
            assert ruszinko_capacity(lower - 1, d) < 1000


def test_ruszinko_capacity_is_monotone_in_rows():
    for d in (2, 3, 5, 10):
        values = [ruszinko_capacity(m, d) for m in range(30)]
        assert values == sorted(values)


def test_weighted_ruszinko_bound_inverts_layer_capacity():
    lower = ruszinko_maximum_weight_lower_bound(1000, 3, 8)
    assert ruszinko_maximum_weight_capacity(lower, 3, 8) >= 1000
    assert ruszinko_maximum_weight_capacity(lower - 1, 3, 8) < 1000
