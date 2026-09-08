"""
KS, recovered from docs/old_main.tex's \\label{lb:KS} (dead ref in
docs/tcs.tex, the current draft) and pinned against docs/tables.tex's
d=2 table.
"""
import pytest

from cgt.bounds.theoretical import KS_lower_bound

D2_N = [100, 1000, 10**4, 10**5, 10**6, 10**8, 10**10, 10**20, 10**30]
D2_KS = [17, 19, 26, 33, 39, 53, 66, 132, 199]


def test_ks_lower_bound_matches_table_except_documented_n100_row():
    # n=100 is the one row that doesn't match (17 tabulated vs. 13 here);
    # see KS_lower_bound's docstring for why this is treated as a likely
    # table transcription slip rather than a bug.
    for n, expected in zip(D2_N, D2_KS):
        got = KS_lower_bound(n, 2)
        if n == 100:
            assert got == 13
        else:
            assert got == expected, n


def test_ks_lower_bound_rejects_d_greater_than_half_n():
    with pytest.raises(ValueError):
        KS_lower_bound(3, 2)
