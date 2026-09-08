from cgt.bounds.balanced_dyachkov_rykov import balanced_partition_capacity
from cgt.bounds.kruskal_katona import kk_capacity, kk_upper_shadow, kk_was_exact


def test_upper_shadow_examples():
    assert kk_upper_shadow(5, 1, 1) == 4
    assert kk_upper_shadow(5, 1, 5) == 10
    assert kk_upper_shadow(5, 2, 10) == 10


def test_kk_refines_balanced_capacity_for_small_parameters():
    assert kk_was_exact(2, 5)
    assert kk_capacity(10, 2, 5) <= float(balanced_partition_capacity(10, 2, 5)) + 1e-9
