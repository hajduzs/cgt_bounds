#!/usr/bin/env python3
"""Generate configurable, transposed construction--bound tables from XML."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import NamedTuple
import xml.etree.ElementTree as ET


DS = (2, 3, 5, 10)
PAIRS = ((2, 3), (5, 10))


class GeneralBound(NamedTuple):
    label: str
    xml_tag: str
    bold: bool = False


class ConstructionQuantity(NamedTuple):
    key: str
    label: str
    xml_suffix: str | None
    bold: bool = False


# Table membership and order are controlled only by these three tuples.
# Adding or removing a bound, construction, or quantity requires no change to
# the rendering logic below.
GENERAL_BOUNDS = (
    GeneralBound(r"$L_{\mathrm F}$", "lower_furedi"),
    GeneralBound(r"$L_{\mathrm R}$", "lower_ruszinko"),
    GeneralBound(r"$L_{\mathrm{DR}}$", "lower_dr"),
    GeneralBound(r"$L_{\mathrm{KS}}$", "lower_ks"),
    GeneralBound(r"$L_{\mathrm{SG}}$", "lower_sg"),
    GeneralBound(r"$\mathbf L_{\infty}$", "lower_infinity", bold=True),
)
METHODS = (
    ("Macula", "macula"), ("EGH-bktrk", "egh_backtrack"),
    ("KS-B", "ks_b"), ("STW-SHF", "stw_shf"), ("KS-C", "ks_c"),
    ("PR", "vpr"), ("GRV", "grv"),
)
BASE_CONSTRUCTION_QUANTITIES = (
    ConstructionQuantity("m", r"$m$", "m"),
    ConstructionQuantity("w", r"$w$", "w"),
    ConstructionQuantity(
        "lbr", r"$L_R$", "lower_ruszinko_weight"
    ),
    ConstructionQuantity(
        "lbdr", r"$L^{\mathrm{direct}}_{\mathrm{DR}}$",
        "lower_dr_weight"
    ),
    ConstructionQuantity("lbsg", r"$L_{\mathrm{SG}}$", "lower_sg_weight"),
    ConstructionQuantity("lbiko", r"$L_{IK\ddot{O}}$", "lower_iko"),
    ConstructionQuantity(
        "lbbal", r"$\mathbf L_{\mathrm{bal}}$", "lower_direct_weight",
        bold=True,
    ),
    ConstructionQuantity(
        "lbrec", r"$\mathbf L_{\mathrm{rec}}$", "lower_recursive", bold=True
    ),
    ConstructionQuantity("gap", r"$\mathrm{GAP}_{\mathrm{best}}$", None),
)


def construction_quantities(include_kk: bool = False) -> tuple[ConstructionQuantity, ...]:
    """Return displayed rows; the expensive KK diagnostic is opt-in."""
    if not include_kk:
        return BASE_CONSTRUCTION_QUANTITIES
    return (
        *BASE_CONSTRUCTION_QUANTITIES[:2],
        ConstructionQuantity("lbkk", r"$L_{KK,w}$", "lower_kk"),
        *BASE_CONSTRUCTION_QUANTITIES[2:],
    )


def applicable_quantities(
    rows: list[ET.Element], prefix: str, *, include_kk: bool = False
) -> tuple[ConstructionQuantity, ...]:
    """Hide an SG row when its hypotheses fail at every displayed point."""
    quantities = construction_quantities(include_kk)
    if any((_integer(row, f"{prefix}_lower_sg_weight") or 0) > 0
           for row in rows):
        return quantities
    return tuple(quantity for quantity in quantities if quantity.key != "lbsg")


def _integer(row: ET.Element, tag: str) -> int | None:
    text = row.findtext(tag)
    return None if text is None else int(text)


def _general_values(row: ET.Element) -> list[int]:
    return [value for bound in GENERAL_BOUNDS
            if (value := _integer(row, bound.xml_tag)) is not None]


def _method_values(
    row: ET.Element, prefix: str, quantities: tuple[ConstructionQuantity, ...]
) -> dict[str, int | float | None]:
    m = _integer(row, f"{prefix}_m")
    if m is None:
        return {quantity.key: None for quantity in quantities}
    values: dict[str, int | float | None] = {"m": m}
    for quantity in quantities:
        if quantity.key in ("m", "gap"):
            continue
        values[quantity.key] = (
            None if quantity.xml_suffix is None
            else _integer(row, f"{prefix}_{quantity.xml_suffix}")
        )
        if quantity.key == "lbsg" and values[quantity.key] == 0:
            values[quantity.key] = None
    applicable = _general_values(row) + [
        int(value) for key, value in values.items()
        if key.startswith("lb") and value is not None
    ]
    values["gap"] = m / max(applicable)
    return values


def _render(value: int | float | None, key: str, bold: bool = False) -> str:
    if value is None:
        rendered = "--"
    elif key == "gap":
        rendered = f"{float(value):.3f}"
    else:
        rendered = str(int(value))
    return rf"\textbf{{{rendered}}}" if bold and value is not None else rendered


def write_table(xml_path: Path, output: Path, *, include_kk: bool = False) -> None:
    rows = ET.parse(xml_path).findall("./measurements/measurement")
    exponents = [int(row.findtext("n_exponent")) for row in rows]
    lines = [
        rf"\begin{{tabular}}{{p{{20mm}}l|{'r' * len(rows)}}}", r"\hline",
        " & ".join(("Construction", "Quantity",
                    *(rf"$10^{{{e}}}$" for e in exponents))) + r" \\",
        r"\hline\hline",
    ]
    for bound_index, bound in enumerate(GENERAL_BOUNDS):
        rendered = [_render(_integer(row, bound.xml_tag), "bound", bound.bold)
                    for row in rows]
        group = (
            rf"\multirow{{{len(GENERAL_BOUNDS)}}}{{20mm}}{{\centering Weight-independent lower bounds}}"
            if bound_index == 0 else ""
        )
        lines.append(" & ".join((group, bound.label, *rendered)) + r" \\")
    lines.append(r"\hline\hline")

    for method_index, (name, prefix) in enumerate(METHODS):
        quantities = applicable_quantities(rows, prefix, include_kk=include_kk)
        values = [_method_values(row, prefix, quantities) for row in rows]
        for quantity_index, quantity in enumerate(quantities):
            rendered = [_render(value[quantity.key], quantity.key, quantity.bold)
                        for value in values]
            lines.append(" & ".join((
                name if quantity_index == 0 else "", quantity.label, *rendered,
            )) + r" \\")
        if method_index != len(METHODS) - 1:
            lines.append(r"\hline")
    lines.extend((r"\hline", r"\end{tabular}"))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_paired_table(
    left_xml: Path,
    right_xml: Path,
    output: Path,
    *,
    left_d: int,
    right_d: int,
    include_kk: bool = False,
) -> None:
    """Put two d-values under adjacent multicolumn headers in one table."""
    left_rows = ET.parse(left_xml).findall("./measurements/measurement")
    right_rows = ET.parse(right_xml).findall("./measurements/measurement")
    left_exponents = [int(row.findtext("n_exponent")) for row in left_rows]
    right_exponents = [int(row.findtext("n_exponent")) for row in right_rows]
    if left_exponents != right_exponents:
        raise ValueError(
            f"paired tables require identical n grids: {left_exponents} != "
            f"{right_exponents}"
        )
    count = len(left_rows)
    lines = [
        rf"\begin{{tabular}}{{p{{20mm}}l|{'r' * count}|{'r' * count}}}",
        r"\hline",
        " & ".join((
            "", "",
            r"\multicolumn{{{}}}{{c|}}{{$d={}$}}".format(count, left_d),
            r"\multicolumn{{{}}}{{c}}{{$d={}$}}".format(count, right_d),
        )) + r" \\",
        " & ".join((
            "Construction", "Quantity",
            *(rf"$10^{{{e}}}$" for e in left_exponents),
            *(rf"$10^{{{e}}}$" for e in right_exponents),
        )) + r" \\",
        r"\hline\hline",
    ]

    for bound_index, bound in enumerate(GENERAL_BOUNDS):
        left = [_render(_integer(row, bound.xml_tag), "bound", bound.bold)
                for row in left_rows]
        right = [_render(_integer(row, bound.xml_tag), "bound", bound.bold)
                 for row in right_rows]
        group = (
            rf"\multirow{{{len(GENERAL_BOUNDS)}}}{{20mm}}{{\centering "
            r"Weight-independent lower bounds}"
            if bound_index == 0 else ""
        )
        lines.append(" & ".join((group, bound.label, *left, *right)) + r" \\")
    lines.append(r"\hline\hline")

    for method_index, (name, prefix) in enumerate(METHODS):
        quantities = applicable_quantities(
            [*left_rows, *right_rows], prefix, include_kk=include_kk
        )
        left_values = [_method_values(row, prefix, quantities)
                       for row in left_rows]
        right_values = [_method_values(row, prefix, quantities)
                        for row in right_rows]
        for quantity_index, quantity in enumerate(quantities):
            left = [_render(value[quantity.key], quantity.key, quantity.bold)
                    for value in left_values]
            right = [_render(value[quantity.key], quantity.key, quantity.bold)
                     for value in right_values]
            lines.append(" & ".join((
                name if quantity_index == 0 else "", quantity.label,
                *left, *right,
            )) + r" \\")
        if method_index != len(METHODS) - 1:
            lines.append(r"\hline")
    lines.extend((r"\hline", r"\end{tabular}"))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--include-kk", action="store_true",
        help="include the optional per-construction Kruskal--Katona row",
    )
    args = parser.parse_args()
    for left_d, right_d in PAIRS:
        write_paired_table(
            args.xml_dir / f"cgt_gap_d{left_d}.xml",
            args.xml_dir / f"cgt_gap_d{right_d}.xml",
            args.output_dir / f"cgt_gap_d{left_d}_d{right_d}.tex",
            left_d=left_d,
            right_d=right_d,
            include_kk=args.include_kk,
        )


if __name__ == "__main__":
    main()
