#!/usr/bin/env python3
"""Generate the article's KS-B lower-bound comparison as PGFPlots LaTeX.

The input is exactly the no-KK measurement XML used by the detailed numerical
tables.  Article notation maps to XML fields as follows:

* ``L_infinity`` is ``lower_infinity``;
* ``L_prev,w`` is the maximum of the previously known weight-aware
  Ruszinko, direct D'yachkov--Rykov, and Inan--Kairouz--Ozgur certificates;
* ``L_w`` is the maximum of ``L_prev,w`` and our recursive weight-aware
  certificate, ``ks_b_lower_recursive``.

The generator checks that the last quantity agrees with the combined KS-B
certificate stored in the XML.  This makes schema or pipeline drift fail
loudly instead of silently changing the plotted curve.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import NamedTuple
import xml.etree.ElementTree as ET


DS = (2, 3, 5, 10)
EXPECTED_EXPONENTS = tuple(range(2, 11))
CONSTRUCTION_OPTIONS = (
    "black,dashed,line width=0.85pt,mark=square*,mark size=1.7pt"
)
INFINITY_OPTIONS = "gray!75!black,solid,line width=0.85pt"
PREVIOUS_OPTIONS = "ksbOrange,solid,line width=0.9pt"
NEW_OPTIONS = "ksbBlue,solid,line width=1.35pt"
WEIGHT_OPTIONS = (
    "ksbGreen,densely dotted,line width=0.9pt,mark=triangle*,mark size=1.9pt"
)
PREVIOUS_WEIGHT_TAGS = (
    "ks_b_lower_ruszinko_weight",
    "ks_b_lower_dr_weight",
    "ks_b_lower_iko",
)


class Point(NamedTuple):
    exponent: int
    m: int
    weight: int
    lower_infinity: int
    lower_previous_weight: int
    lower_weight: int


def _required_int(row: ET.Element, tag: str, source: Path) -> int:
    text = row.findtext(tag)
    if text is None:
        raise ValueError(f"{source}: missing required XML field <{tag}>")
    return int(text)


def read_ksb_points(xml_dir: Path, d: int) -> list[Point]:
    """Read and validate one no-KK KS-B measurement series."""
    source = xml_dir / f"cgt_gap_d{d}.xml"
    rows = ET.parse(source).findall("./measurements/measurement")
    if any(row.find("ks_b_lower_kk") is not None for row in rows):
        raise ValueError(
            f"{source}: contains KK measurements; regenerate without --include-kk"
        )

    points = []
    for row in rows:
        exponent = _required_int(row, "n_exponent", source)
        previous = max(_required_int(row, tag, source)
                       for tag in PREVIOUS_WEIGHT_TAGS)
        recursive = _required_int(row, "ks_b_lower_recursive", source)
        lower_weight = max(previous, recursive)
        stored = _required_int(row, "ks_b_lower_combined_weight", source)
        if lower_weight != stored:
            raise ValueError(
                f"{source}, n=10^{exponent}: computed L_w={lower_weight}, "
                f"but the detailed-table field stores {stored}"
            )
        point = Point(
            exponent=exponent,
            m=_required_int(row, "ks_b_m", source),
            weight=_required_int(row, "ks_b_w", source),
            lower_infinity=_required_int(row, "lower_infinity", source),
            lower_previous_weight=previous,
            lower_weight=lower_weight,
        )
        if not (
            point.lower_previous_weight <= point.lower_weight <= point.m
            and point.lower_infinity <= point.m
        ):
            raise ValueError(f"{source}, n=10^{exponent}: inconsistent KS-B bounds")
        points.append(point)

    exponents = tuple(point.exponent for point in points)
    if exponents != EXPECTED_EXPONENTS:
        raise ValueError(
            f"{source}: expected exponents {EXPECTED_EXPONENTS}, got {exponents}"
        )
    return points


def _coordinates(points: list[Point], field: str) -> str:
    return " ".join(
        f"({point.exponent},{getattr(point, field)})" for point in points
    )


def _primary_panel(points: list[Point], d: int, panel: int) -> list[str]:
    ylabel = (
        r"ylabel={Number of tests $m$},"
        if panel in (1, 3) else ""
    )
    legend = [
        "legend to name={ksb-weight-aware-legend},",
        r"\addlegendentry{KS-B construction $m$}",
        r"\addlegendentry{$L_\infty$}",
        r"\addlegendentry{$L_{\mathrm{ref},w}$}",
        r"\addlegendentry{$L_w$}",
        rf"\addlegendimage{{{WEIGHT_OPTIONS}}}",
        r"\addlegendentry{Column weight $w$}",
    ] if panel == 1 else []
    return [
        rf"\nextgroupplot[title={{$d={d}$}},{ylabel}",
        *(legend[:1]),
        "]",
        rf"\addplot[{CONSTRUCTION_OPTIONS}] coordinates {{{_coordinates(points, 'm')}}};",
        *(legend[1:2]),
        rf"\addplot[{INFINITY_OPTIONS}] coordinates {{{_coordinates(points, 'lower_infinity')}}};",
        *(legend[2:3]),
        rf"\addplot[{PREVIOUS_OPTIONS}] coordinates {{{_coordinates(points, 'lower_previous_weight')}}};",
        *(legend[3:4]),
        rf"\addplot[{NEW_OPTIONS}] coordinates {{{_coordinates(points, 'lower_weight')}}};",
        *(legend[4:]),
    ]


def render_figure(series: dict[int, list[Point]]) -> str:
    """Render one compact 2x2 PGFPlots figure body."""
    lines = [
        r"% Generated by scripts/generate_ksb_weight_aware_bounds_figure.py.",
        r"% Input: output/article-measurements/cgt_gap_d{2,3,5,10}.xml (no KK).",
        r"% TikZ/PGFPlots fragment; include it with \input from the article.",
        r"\definecolor{ksbBlue}{RGB}{0,92,171}",
        r"\definecolor{ksbOrange}{RGB}{230,126,34}",
        r"\definecolor{ksbGreen}{RGB}{0,128,96}",
        r"\begin{tikzpicture}",
        r"\begin{groupplot}[",
        r"  group style={group size=2 by 2,horizontal sep=11mm,vertical sep=12mm},",
        r"  width=0.375\textwidth,height=0.25\textwidth,scale only axis,",
        r"  xmin=1,xmax=10,xtick={1,...,10},",
        r"  xticklabels={$10^1$,$10^2$,$10^3$,$10^4$,$10^5$,$10^6$,$10^7$,$10^8$,$10^9$,$10^{10}$},",
        r"  x tick label style={font=\scriptsize,rotate=45,anchor=east},",
        r"  ymode=log,log basis y=10,",
        r"  tick label style={font=\scriptsize},label style={font=\small},",
        r"  title style={font=\small},grid=major,grid style={gray!18},",
        r"  legend style={font=\scriptsize,draw=none,legend columns=5,/tikz/every even column/.append style={column sep=4pt}},",
        r"]",
    ]
    for panel, d in enumerate(DS, start=1):
        lines.extend(_primary_panel(series[d], d, panel))
    lines.append(r"\end{groupplot}")

    for panel, d in enumerate(DS, start=1):
        column = 1 if panel in (1, 3) else 2
        row = 1 if panel in (1, 2) else 2
        ylabel = (
            r"ylabel={Column weight $w$}," if column == 2 else r"ylabel={},"
        )
        lines.extend((
            r"\begin{axis}[",
            rf"  at={{(group c{column}r{row}.south west)}},anchor=south west,",
            r"  width=0.375\textwidth,height=0.25\textwidth,scale only axis,",
            r"  xmin=1,xmax=10,axis x line=none,axis y line*=right,",
            rf"  {ylabel}",
            r"  ymin=0,ymax=50,ytick={0,10,...,50},",
            r"  tick label style={font=\scriptsize},label style={font=\small},",
            r"  grid=none,",
            r"]",
            rf"\addplot[{WEIGHT_OPTIONS}] coordinates {{{_coordinates(series[d], 'weight')}}};",
            r"\end{axis}",
        ))
    lines.extend((
        r"\node[anchor=north] at ([yshift=-4mm]current bounding box.south) "
        r"{\ref{ksb-weight-aware-legend}};",
        r"\end{tikzpicture}",
    ))
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    series = {d: read_ksb_points(args.xml_dir, d) for d in DS}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_figure(series), encoding="utf-8")


if __name__ == "__main__":
    main()
