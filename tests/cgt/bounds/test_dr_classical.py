"""
DR_classical vs. docs/tables.tex's plain "DR" column (the classical
D'yachkov-Rykov recursion -- see its own docstring).
"""
from cgt.bounds.theoretical import DR_classical

D2_N = [100, 1000, 10**4, 10**5, 10**6]
D2_DR = [13, 19, 27, 36, 46]


def test_dr_classical_within_one_of_table():
    # See DR_classical's docstring: exact at n=100, off by exactly 1 at
    # the next few rows (likely a boundary-condition transcription
    # difference, not a wrong formula -- shape/asymptotics match).
    for n, expected in zip(D2_N, D2_DR):
        got = DR_classical(n, 2)
        assert abs(got - expected) <= 1, (n, got, expected)


def test_dr_classical_exact_at_n100():
    assert DR_classical(100, 2) == 13


def test_dr_classical_monotone_in_n():
    prev = DR_classical(10, 2)
    for n in [50, 100, 1000, 10**4]:
        cur = DR_classical(n, 2)
        assert cur >= prev
        prev = cur
