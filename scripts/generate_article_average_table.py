#!/usr/bin/env python3
"""Transpose the article's construction summary and report all three LBs."""

from __future__ import annotations

import argparse
from math import ceil, log
from pathlib import Path
import statistics
import warnings
import xml.etree.ElementTree as ET


DS = (2, 3, 5, 10)
METHODS = (
    "KS-A/CWC", "KS-B", "STW-SHF", "KS-C", "HS", "Macula", "GRV",
    "EGH-gen", "EGH-bktrk", "PR",
)
PREFIXES = {
    "KS-B": "ks_b", "STW-SHF": "stw_shf", "KS-C": "ks_c",
    "Macula": "macula", "GRV": "grv", "EGH-gen": "egh_general",
    "EGH-bktrk": "egh_backtrack", "PR": "vpr",
}


def _mean(items: list[float | None]) -> float | None:
    values = [item for item in items if item is not None]
    return statistics.fmean(values) if values else None


def _get(row: ET.Element, tag: str, cast, missing: set[str]) -> object | None:
    """Read ``tag`` from ``row``, casting with ``cast``.

    Measurement XML is regenerated incrementally, so it's normal for older
    files to be missing fields a newer script version added. Rather than
    crash, record the gap in ``missing`` (surfaced once, as a warning, by
    the caller) and let downstream code treat the value as unavailable.
    """
    text = row.findtext(tag)
    if text is None:
        missing.add(tag)
        return None
    return cast(text)


def _general(row: ET.Element, missing: set[str]) -> int | None:
    values = [
        v for v in (
            _get(row, tag, int, missing) for tag in
             ("lower_furedi", "lower_ruszinko", "lower_dr", "lower_ks",
              "lower_sg", "lower_infinity")
        ) if v is not None
    ]
    return max(values) if values else None


def averages(xml_dir: Path) -> dict[str, dict[int, dict[str, float | None]]]:
    """Return mean LBs and the mean per-instance best-LB gap.

    Any measurement field missing from the input XML (e.g. an incompletely
    regenerated table) is treated as unavailable rather than fatal: the
    affected mean/gap entries fall back to ``None`` (rendered as ``--``),
    and a single warning at the end lists every field that was missing
    somewhere, so the incompleteness is visible without aborting the run.
    """
    missing: set[str] = set()
    result = {name: {} for name in METHODS}
    for d in DS:
        rows = ET.parse(xml_dir / f"cgt_gap_d{d}.xml").findall(
            "./measurements/measurement"
        )
        for name, prefix in PREFIXES.items():
            samples: list[tuple[int, int | None, int | None, int | None, float | None]] = []
            for row in rows:
                m_text = row.findtext(f"{prefix}_m")
                if m_text is None:
                    continue
                lb0 = _general(row, missing)
                # ``lower_weight`` is the backward-compatible fallback for
                # XML files generated before the raw maximum-weight field was
                # introduced.  Fresh measurements always use the first tag;
                # only flag it missing if *both* are absent.
                lbmax_text = (row.findtext(f"{prefix}_lower_maximum")
                              or row.findtext(f"{prefix}_lower_weight"))
                if lbmax_text is None:
                    missing.add(f"{prefix}_lower_maximum")
                lbmax = int(lbmax_text) if lbmax_text is not None else None
                lbavg = _get(row, f"{prefix}_lower_average", int, missing)
                # Recompute the ratio from exactly the certificate rows shown
                # in the detailed tables.  Do not trust a cached XML ratio:
                # that could silently make the summary disagree with Tables
                # 3--4 after the displayed certificate set changes.
                displayed_weight_bounds = [
                    _get(row, f"{prefix}_{suffix}", int, missing)
                    for suffix in (
                        "lower_ruszinko_weight", "lower_dr_weight",
                        "lower_sg_weight", "lower_iko", "lower_direct_weight",
                        "lower_recursive",
                    )
                ]
                certificates = [
                    value for value in (lb0, *displayed_weight_bounds)
                    if value is not None
                ]
                m = int(m_text)
                gap = m / max(certificates) if certificates else None
                samples.append((m, lb0, lbmax, lbavg, gap))
            result[name][d] = {
                "lb0": _mean([x[1] for x in samples]),
                "lbmax": _mean([x[2] for x in samples]),
                "lbavg": _mean([x[3] for x in samples]),
                # The XML GAP uses the strongest certificate reported in the
                # detailed table. Average ratios, never averaged numerators
                # and denominators separately.
                "gap": _mean([x[4] for x in samples]),
            }

        general = [_general(row, missing) for row in rows]
        hs_m = []
        for row in rows:
            n = 10 ** int(row.findtext("n_exponent"))
            hs_m.append(ceil(16 * d * d *
                             (1 + log(2, 3) + log(2, 3) * log(n, 2))))
        result["HS"][d] = {
            "lb0": _mean(general), "lbmax": None, "lbavg": None,
            "gap": _mean([
                m / lb if lb is not None else None
                for m, lb in zip(hs_m, general, strict=True)
            ]),
        }
    cwc_rows = ET.parse(xml_dir / "cwc_small.xml").findall(
        "./measurements/measurement"
    )
    for d in DS:
        samples = [row for row in cwc_rows if int(row.findtext("d")) == d]
        lb0s = [_get(row, "lower_general", int, missing) for row in samples]
        lbmaxs = []
        for row in samples:
            text = row.findtext("lower_maximum") or row.findtext("lower_weight")
            if text is None:
                missing.add("lower_maximum")
            lbmaxs.append(int(text) if text is not None else None)
        lbavgs = [_get(row, "lower_average", int, missing) for row in samples]
        gaps = []
        for row, lb0, lbmax, lbavg in zip(samples, lb0s, lbmaxs, lbavgs, strict=True):
            candidates = [v for v in (lb0, lbmax, lbavg) if v is not None]
            gaps.append(int(row.findtext("m")) / max(candidates) if candidates else None)
        result["KS-A/CWC"][d] = {
            "lb0": _mean(lb0s), "lbmax": _mean(lbmaxs), "lbavg": _mean(lbavgs),
            "gap": _mean(gaps),
        }

    if missing:
        warnings.warn(
            "Incomplete article-measurement XML in "
            f"{xml_dir}: missing field(s) {sorted(missing)} in one or more "
            "rows. Affected table entries were rendered as '--' instead of "
            "aborting -- regenerate the measurements "
            "(scripts/generate_article_measurements.sh) for a complete table.",
            stacklevel=2,
        )
    return result


def _render(value: float | None, *, gap: bool = False) -> str:
    if value is None:
        return "--"
    if not gap:
        return f"{value:.1f}"
    magnitude = abs(value)
    if magnitude < 10:
        return f"{value:.3f}"
    if magnitude < 100:
        return f"{value:.2f}"
    if magnitude < 1000:
        return f"{value:.1f}"
    return f"{value:.0f}"


def write_latex(
    summary: dict[str, dict[int, dict[str, float | None]]], output: Path
) -> None:
    """Write the literature-style table containing only mean best gaps."""
    metadata = {
"KS-A/CWC": (r"\cite[Sec.~V-A,D]{kautz1964nonrandom}",
             "1964", "", r"Best BCW-ECC from Brouwer's tables."),
"KS-B": (r"\cite[Sec.~V-B]{kautz1964nonrandom}", "1964", r"\cmark",
         "RS code with one-hot expansion."),
"STW-SHF": (r"\cite[Thm.~4.8]{stinson2000secure}", "2000", r"\cmark",
            "Recursive separating hash family."),
"KS-C": (r"\cite[Sec.~V-C]{kautz1964nonrandom}", "1964", r"\cmark",
         "Recursive ECC construction."),
"HS": (r"\cite{hwang1987non}", "1987", "",
       "Greedy upper estimate reported by EGH."),
"Macula": (r"\cite{macula1996simple}", "1996", "",
           "Subset--inclusion construction."),
"GRV": (r"\cite{gargano2020low}", "2020", "",
        r"LLL constant-weight construction."),
"EGH-gen": (r"\cite[Sec.~2]{eppstein2007improved}", "2007", r"\cmark",
            "General CRT sieve."),
"EGH-bktrk": (r"\cite[Sec.~2]{eppstein2007improved}", "2007", r"\cmark",
              "Backtracking CRT sieve."),
"PR": (r"\cite{porat2011explicit}", "2011", r"\cmark",
       "ECC meeting the GV bound."),
    }
    lines = [
        r"\begin{tabular}{c|c|c|c|rrrr|l}", r"\hline",
        r"Abbrev. & Citation & Year & UB & \multicolumn{4}{c|}{mean $\mathrm{GAP}_{\mathrm{best}}$} & Remark \\",
        r" & & & & $d=2$ & $d=3$ & $d=5$ & $d=10$ & \\",
        r"\hline\hline",
    ]
    for name, fields in metadata.items():
        rendered = [_render(summary[name][d]["gap"], gap=True) for d in DS]
        lines.append(" & ".join((name, *fields[:3], *rendered, fields[3])) + r" \\")
    lines.extend((r"\hline", r"\end{tabular}"))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_latex(averages(args.xml_dir), args.output)


if __name__ == "__main__":
    main()
