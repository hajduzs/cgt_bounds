"""
Pins cgt.bounds.theoretical's new balanced-partition certificates against
the numbers actually printed in docs/tcs.tex (see M_bal_unrestricted's and
M_BDR's docstrings for the exact provenance of each).
"""
from cgt.bounds.theoretical import U_bal, M_bal_unrestricted, M_BDR


def test_u_bal_matches_d2_companion_table_LB_star():
    expected = {
        100: 13, 1000: 22, 10**4: 32, 10**5: 42, 10**6: 52,
        10**8: 71, 10**10: 92, 10**20: 193, 10**30: 295,
    }
    for n, m in expected.items():
        assert M_bal_unrestricted(n, 2) == m, n


def test_u_bal_monotone_in_m():
    for m in range(1, 40):
        assert U_bal(m, 2, m) <= U_bal(m + 1, 2, m + 1)


def test_m_bdr_dominates_sperner_base_case():
    # M_BDR(n, 1) must equal the Sperner bound used as its own base case.
    import math

    def sperner_m(n):
        m = 1
        while math.comb(m, m // 2) < n:
            m += 1
        return m

    for n in [10, 100, 1000]:
        assert M_BDR(n, 1) == sperner_m(n)
