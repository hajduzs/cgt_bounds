#!/usr/bin/env python3
"""Generate the small-length constant-weight-code evaluation XML."""

from __future__ import annotations

import argparse
from math import floor
from pathlib import Path
import xml.etree.ElementTree as ET

from cgt.bounds.balanced_dyachkov_rykov import (
    average_weight_profile_relaxation_lower_bound,
    combined_average_weight_lower_bound,
    balanced_dr_lower_bound,
    recursive_maximum_weight_lower_bound,
)
from cgt.bounds.classical import dyachkov_rykov_lower_bound, kautz_singleton_lower_bound
from cgt.bounds.furedi import furedi_lower_bound
from cgt.bounds.theoretical import D_star
from cgt.bounds.kruskal_katona import kk_maximum_weight_lower_bound, kk_was_exact
from cgt.bounds.literature import (
    inan_kairouz_ozgur_lower_bound,
    shangguan_ge_lower_bound,
    shangguan_ge_weight_aware_lower_bound,
)
from cgt.bounds.ruszinko import ruszinko_lower_bound
from cwc.bounds.brouwer import loadBounds


DEFECTIVES = (2, 3, 5)
LENGTHS = (10, 15, 20, 25, 30, 40)


def _text(parent: ET.Element, name: str, value: object) -> None:
    ET.SubElement(parent, name).text = str(value)


def best_cwc(bounds, m: int, d: int):
    """Return the largest stored Brouwer lower bound satisfying D >= D*(d,w)."""
    candidates = []
    for (length, distance, weight), (lower, _upper, description) in bounds.items():
        if length == m and lower is not None and distance >= D_star(d, weight):
            candidates.append((int(lower), -weight, -distance, weight, distance, description))
    return max(candidates) if candidates else None


def _johnson_cwc_upper_bound(m: int, distance: int, weight: int) -> int:
    """The iterated Johnson upper bound for ``A(m,distance,weight)``."""
    value = 1
    for overlap in range(weight - distance // 2, -1, -1):
        value = floor(value * (m - overlap) / (weight - overlap))
    return value


def cwc_table_row_lower_bound(bounds, n: int, d: int, weight: int) -> int:
    """Invert Johnson plus applicable stored Brouwer/AVZ CWC upper bounds.

    A stored upper bound at distance ``D <= D_star(d, weight)`` is also a
    valid upper bound at the required (larger) distance.  The scan stops at
    the first length not ruled out, so missing table entries cannot create an
    unjustified jump.
    """
    required_distance = D_star(d, weight)
    m = weight
    while True:
        upper = _johnson_cwc_upper_bound(m, required_distance, weight)
        stored = [
            int(candidate_upper)
            for (length, distance, candidate_weight),
                (_lower, candidate_upper, _description) in bounds.items()
            if length == m and candidate_weight == weight
            and distance <= required_distance and candidate_upper is not None
        ]
        if stored:
            upper = min(upper, *stored)
        if upper >= n:
            return m
        m += 1


def _applicable_general_lower_bound(n: int, d: int) -> int:
    values = [balanced_dr_lower_bound(n, d), shangguan_ge_lower_bound(n, d)]
    if d >= 2:
        values.append(ruszinko_lower_bound(n, d))
    for bound in (
        furedi_lower_bound,
        dyachkov_rykov_lower_bound,
        kautz_singleton_lower_bound,
    ):
        try:
            values.append(bound(n, d))
        except ValueError:
            pass
    return max(values)


def measurement_xml(*, include_kk: bool = False) -> ET.Element:
    bounds, _ = loadBounds()
    root = ET.Element("results")
    _text(root, "version", "1.0")
    args = ET.SubElement(root, "args")
    _text(args, "par_command", "article-small-cwc-measurement")
    measurements = ET.SubElement(root, "measurements")
    row_id = 0
    for d in DEFECTIVES:
        for m in LENGTHS:
            best = best_cwc(bounds, m, d)
            if best is None:
                continue
            n, _neg_w, _neg_D, weight, distance, description = best
            recursive_lower = recursive_maximum_weight_lower_bound(n, d, weight)
            kk_lower = kk_maximum_weight_lower_bound(n, d, weight) if include_kk else 0
            iko_lower = inan_kairouz_ozgur_lower_bound(n, d, weight)
            sg_weight_lower = shangguan_ge_weight_aware_lower_bound(n, d, weight)
            cwc_table_lower = cwc_table_row_lower_bound(bounds, n, d, weight)
            lower = max(
                recursive_lower, kk_lower, iko_lower, sg_weight_lower,
                cwc_table_lower,
            )
            lower_general = _applicable_general_lower_bound(n, d)
            lower_average = combined_average_weight_lower_bound(
                n, d, weight, lower_general
            )
            lower_average_profile = max(
                lower_general,
                average_weight_profile_relaxation_lower_bound(n, d, weight),
            )
            lower_maximum = max(lower, lower_general)
            lower = max(lower_maximum, lower_average)
            gap = m / lower
            # Table 4 is intended as a compact list of near-optimal small
            # CWC constructions; less informative rows are omitted.
            if gap > 2 + 1e-12:
                continue
            row = ET.SubElement(measurements, "measurement")
            row_id += 1
            for name, value in (
                ("row_id", row_id),
                ("d", d), ("m", m), ("weight", weight),
                ("distance", distance), ("size", n),
                ("lower_maximum", lower_maximum),
                ("lower_weight", lower), ("gap", f"{gap:.6f}"),
                ("lower_average", lower_average),
                ("lower_average_profile", lower_average_profile),
                ("gap_average", f"{m / lower_average:.6f}"),
                ("lower_general", lower_general),
                ("gap_unweighted", f"{m / lower_general:.6f}"),
                ("lower_iko", iko_lower),
                ("lower_sg_weight", sg_weight_lower),
                ("lower_cwc_table", cwc_table_lower),
                ("source", description),
            ):
                _text(row, name, value)
            if include_kk:
                _text(row, "lower_kk", kk_lower)
                _text(row, "kk_exact", int(kk_was_exact(d, weight)))
    return root


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--include-kk", action="store_true",
        help="include Kruskal--Katona bounds in the best-bound calculation",
    )
    args = parser.parse_args()
    root = measurement_xml(include_kk=args.include_kk)
    ET.indent(root, space="  ")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(args.output, encoding="utf-8", xml_declaration=True)


if __name__ == "__main__":
    main()
