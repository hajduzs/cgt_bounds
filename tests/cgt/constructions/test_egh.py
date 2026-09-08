from cgt.constructions.egh import egh_crt_zfd, egh_parameters
from cgt.constructions.ks64 import is_d_disjunct


def test_egh_paper_table_reproduction():
    assert egh_parameters(100, 5, "general").m == 160
    assert egh_parameters(100, 5, "backtrack").m == 131
    assert egh_parameters(100, 10, "general").m == 440
    assert egh_parameters(100, 10, "backtrack").m == 378


def test_egh_survey_d2_rows():
    ns = (100, 10**3, 10**4, 10**5, 10**6, 10**8, 10**10, 10**20, 10**30)
    assert [egh_parameters(n, 2, "general").m for n in ns] == [
        41, 77, 100, 160, 197, 281, 440, 1264, 2584
    ]
    assert [egh_parameters(n, 2, "backtrack").m for n in ns] == [
        36, 60, 89, 131, 168, 268, 378, 1176, 2350
    ]


def test_materialized_egh_code_is_disjunct_and_has_expected_distance():
    code = egh_crt_zfd(30, 2, "backtrack")
    assert is_d_disjunct(code)
    distances = [
        len(left ^ right)
        for index, left in enumerate(code.columns)
        for right in code.columns[index + 1 :]
    ]
    assert min(distances) == code.metadata["minimum_distance"]
