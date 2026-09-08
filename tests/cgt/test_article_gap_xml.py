import importlib.util
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest


def _load_module():
    script = Path(__file__).parents[2] / "scripts" / "generate_article_gap_xml.py"
    spec = importlib.util.spec_from_file_location("generate_article_gap_xml", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _assert_row_invariants(row, module):
    for tag in ("lower_ruszinko", "lower_dr", "lower_ks"):
        assert int(row.findtext(tag)) > 0
    for prefix in ("egh_general", "egh_backtrack", "macula", "ks_b", "stw_shf", "ks_c", "vpr", "grv"):
        assert float(row.findtext(f"{prefix}_gap")) >= 1
        general = max(int(row.findtext(tag)) for tag in
                      ("lower_furedi", "lower_ruszinko", "lower_dr",
                       "lower_ks", "lower_sg", "lower_infinity"))
        maximum = int(row.findtext(f"{prefix}_lower_maximum"))
        for suffix in ("lower_recursive", "lower_dr_weight",
                       "lower_ruszinko_weight", "lower_direct_weight",
                       "lower_combined_weight", "lower_iko"):
            assert int(row.findtext(f"{prefix}_{suffix}")) > 0
        # lower_sg_weight is 0 (a no-op) whenever neither Shangguan--Ge
        # weight-aware term applies at this construction's weight, so it is
        # not asserted positive like the other weight-aware certificates.
        assert int(row.findtext(f"{prefix}_lower_sg_weight")) >= 0
        # lower_combined_weight is the max over only the weight-aware
        # certificate family (infinity/kk/recursive/direct/iko/sg_weight) --
        # it is not required to dominate `general` (the max over the
        # classical family). The two families are combined by taking the
        # max over every individual certificate (see `expected` below), not
        # by requiring one family's aggregate to dominate the other's.
        assert maximum >= general
        expected = int(row.findtext(f"{prefix}_m")) / max(
            general,
            int(row.findtext(f"{prefix}_lower_recursive")),
            int(row.findtext(f"{prefix}_lower_ruszinko_weight")),
            int(row.findtext(f"{prefix}_lower_dr_weight")),
            int(row.findtext(f"{prefix}_lower_direct_weight")),
            int(row.findtext(f"{prefix}_lower_iko")),
            int(row.findtext(f"{prefix}_lower_sg_weight")),
        )
        assert abs(float(row.findtext(f"{prefix}_gap")) - expected) < 1e-6
        assert row.find(f"{prefix}_lower_kk") is None


def test_generated_gap_xml_fast(tmp_path, monkeypatch):
    """Smoke-checks the same self-consistency invariants as the full sweep,
    but over two small exponents for one d value instead of the full
    DEFECTIVES x TARGET_EXPONENTS grid. The invariants exercised here don't
    depend on which (d, n) cell is checked -- they either hold for the
    formula in general or they don't -- so this is a fast substitute for
    local iteration; see test_generated_gap_xml_full_sweep (slow) for the
    complete grid.

    Two things dominate cost here, so both are kept small: exponent (several
    bound formulas cost ~10x more per step, so n=10**10 alone would swamp a
    "fast" test) and d (kk_upper_shadow/kk_profiles in
    cgt/bounds/kruskal_katona.py do real combinatorial profile enumeration
    that scales steeply with d -- d=10 alone costs ~8x what d=2 does at the
    same exponents). Sticking to d=2 (the cheapest DEFECTIVES value) and
    exponents (2, 5) keeps this under ~5s.
    """
    module = _load_module()
    monkeypatch.setattr(module, "TARGET_EXPONENTS", (2, 5))
    monkeypatch.setattr(module, "TARGETS", tuple(10**e for e in (2, 5)))

    output = tmp_path / "d2.xml"
    module.write_xml(module.measurement_xml(2), output)
    rows = ET.parse(output).findall("./measurements/measurement")
    assert len(rows) == 2
    assert [int(row.findtext("n_exponent")) for row in rows] == [2, 5]
    for row in rows:
        _assert_row_invariants(row, module)


@pytest.mark.slow
def test_generated_gap_xml_full_sweep(tmp_path):
    """Full DEFECTIVES x TARGET_EXPONENTS grid -- the real, complete article
    gap-table generation. Slow (~4.5 min); run explicitly (`pytest -m slow`)
    or in CI, not by default. See test_generated_gap_xml_fast for the
    every-day version of these same checks.
    """
    module = _load_module()

    for d in module.DEFECTIVES:
        output = tmp_path / f"d{d}.xml"
        module.write_xml(module.measurement_xml(d), output)
        rows = ET.parse(output).findall("./measurements/measurement")
        assert len(rows) == len(module.TARGETS)
        assert [int(row.findtext("n_exponent")) for row in rows] == list(range(2, 11))
        for row in rows:
            _assert_row_invariants(row, module)
