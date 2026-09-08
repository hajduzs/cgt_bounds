import warnings
import pytest
from cwc.bounds.brouwer import boundA, boundADesc, boundblockA

@pytest.mark.bounds
def test_bounds_brouwer_optimal():
    """
    Certified-optimal entries (lb == ub). These don't drift when
    data/bounds/brouwer_main.txt is refreshed -- "optimal" doesn't get
    un-proven -- so it's safe to hard-assert the exact numeric bound.
    """
    assert boundA(10, 4, 3) == (13, 13)
    assert boundA(10, 4, 4) == (30, 30)
    assert boundA(10, 4, 5) == (36, 36)
    assert boundA(10, 6, 4) == (5, 5)
    assert boundA(10, 6, 5) == (6, 6)
    assert boundA(20, 10, 7) == (10, 10)
    assert boundA(11, 2, 2) == (55, 55)
    assert boundA(77, 16, 9) == (None, None)

@pytest.mark.bounds
def test_bounds_brouwer_open_problems_sanity():
    """
    Still-open (lb != ub) entries: as better constructions/bounds are
    found, these numbers move. Hard-coding an exact snapshot here would
    make the suite fail every time Brouwer's table is refreshed for
    reasons that have nothing to do with a regression in this codebase --
    only assert the invariant that must always hold.
    """
    for m, D, w in [(23, 6, 10), (19, 6, 5), (39, 6, 6), (56, 8, 7), (64, 12, 7), (69, 14, 8)]:
        lb, ub = boundA(m, D, w)
        if lb is not None and ub is not None:
            assert 0 < lb <= ub, f"A({m},{D},{w}): lb={lb}, ub={ub}"

@pytest.mark.bounds
def test_bounds_brouwer_description_flags_drift():
    """
    boundADesc() returns Brouwer's free-text annotation (e.g.
    "lost:2970-7521" for provenance-unclear entries). That text changes
    whenever the local brouwer_main.txt snapshot is refreshed from the
    live page -- so on a mismatch, warn instead of failing: it signals
    the table moved, not that this codebase broke.
    """
    expected_bound = (2969, 7521)
    assert boundA(23, 6, 10) == expected_bound

    expected_desc = 'lost:2970-7521'
    actual_desc = boundADesc(23, 6, 10)
    if actual_desc != expected_desc:
        warnings.warn(
            f"boundADesc(23,6,10) is {actual_desc!r}, expected {expected_desc!r} -- "
            "data/bounds/brouwer_main.txt was likely refreshed from Brouwer's live "
            "page since this snapshot was taken. This is expected drift, not "
            "necessarily a regression, but worth a quick look if unexpected."
        )

@pytest.mark.bounds
def test_brouwer_block_bounds():
    assert boundblockA(28, 2, 12) == (178, 288)
    assert boundblockA(16, 3, 8) == (387, 1923)
    assert boundblockA(16, 3, 15) == (3, 3)
    assert boundblockA(12, 4, 11) == (4, 4)
    assert boundblockA(11, 5, 9) == (25, 35)
    assert boundblockA(5, 5, 3) == (125, 125)
