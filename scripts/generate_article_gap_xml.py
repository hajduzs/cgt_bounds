#!/usr/bin/env python3
"""Generate the survey's reproducible construction-versus-bound XML tables."""

from __future__ import annotations

import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

from cgt.bounds.balanced_dyachkov_rykov import (
    average_weight_profile_relaxation_lower_bound,
    combined_average_weight_lower_bound,
    balanced_dr_lower_bound,
    maximum_weight_lower_bound,
    recursive_maximum_weight_lower_bound,
)
from cgt.bounds.furedi import furedi_lower_bound
from cgt.bounds.kruskal_katona import (
    kk_maximum_weight_lower_bound,
    kk_unrestricted_lower_bound,
    kk_was_exact,
)
from cgt.bounds.ruszinko import (
    ruszinko_lower_bound,
    ruszinko_maximum_weight_lower_bound,
)
from cgt.bounds.classical import (
    dyachkov_rykov_lower_bound,
    kautz_singleton_lower_bound,
    maximum_weight_dyachkov_rykov_lower_bound,
)
from cgt.bounds.literature import (
    inan_kairouz_ozgur_lower_bound,
    shangguan_ge_lower_bound,
    shangguan_ge_weight_aware_lower_bound,
)
from cgt.constructions.egh import egh_parameters
from cgt.constructions.gargano import gargano_parameters
from cgt.constructions.macula import macula_parameters
from cgt.constructions.stw_shf import stw_parameters
from cgt.tables.inr import inr_parameters
from cgt.tables.ks import _is_prime_power, ks_b_parameters, ks_c_parameters
from cgt.tables.vpr import vpr_parameters


TARGET_EXPONENTS = tuple(range(2, 11))
TARGETS = tuple(10**exponent for exponent in TARGET_EXPONENTS)
DEFECTIVES = (2, 3, 5, 10)


def _text(parent: ET.Element, name: str, value: object) -> None:
    ET.SubElement(parent, name).text = str(value)


def _best_vpr(n: int, d: int) -> tuple[int, int]:
    candidates = []
    for q in range(2, 258):
        if not _is_prime_power(q):
            continue
        try:
            m, weight = vpr_parameters(n, d, q)
        except ValueError:
            continue
        candidates.append((m, weight))
    if not candidates:
        raise RuntimeError(f"no VPR parameter found for n={n}, d={d}")
    return min(candidates)


def construction_rows(n: int, d: int) -> list[tuple[str, int, int]]:
    rows = []
    for variant, label in (("general", "EGH-general"), ("backtrack", "EGH-backtrack")):
        parameters = egh_parameters(n, d, variant)
        rows.append((label, parameters.m, parameters.w))
    macula_m, macula_w, _v, _k = macula_parameters(n, d)
    rows.append(("Macula", macula_m, macula_w))
    ks_m, ks_q = ks_b_parameters(n, d)
    rows.append(("KS-B", ks_m, ks_m // ks_q))
    stw_m, stw_w, _stw_q, _stw_k, _stw_j, _stw_capacity = stw_parameters(n, d)
    rows.append(("STW-SHF", stw_m, stw_w))
    ks_c_m, ks_c_w = ks_c_parameters(n, d)
    rows.append(("KS-C", ks_c_m, ks_c_w))
    vpr_m, vpr_w = _best_vpr(n, d)
    rows.append(("VPR", vpr_m, vpr_w))
    gargano_m, gargano_w = gargano_parameters(n, d)
    rows.append(("GRV", gargano_m, gargano_w))
    return rows


def measurement_xml(d: int, *, include_kk: bool = False) -> ET.Element:
    root = ET.Element("results")
    _text(root, "version", "1.0")
    args = ET.SubElement(root, "args")
    _text(args, "par_command", "article-gap-measurement")
    _text(args, "par_d", d)
    measurements = ET.SubElement(root, "measurements")
    for exponent, n in zip(TARGET_EXPONENTS, TARGETS, strict=True):
        furedi = furedi_lower_bound(n, d)
        ruszinko = ruszinko_lower_bound(n, d)
        dr = dyachkov_rykov_lower_bound(n, d)
        ks = kautz_singleton_lower_bound(n, d)
        sg = shangguan_ge_lower_bound(n, d)
        bdr = balanced_dr_lower_bound(n, d)
        kk_infinity = kk_unrestricted_lower_bound(n, d) if include_kk else 0
        infinity = max(bdr, kk_infinity)
        general = max(furedi, ruszinko, dr, ks, sg, infinity)
        row = ET.SubElement(measurements, "measurement")
        _text(row, "n_exponent", exponent)
        for construction, m, weight in construction_rows(n, d):
            recursive_bound = recursive_maximum_weight_lower_bound(n, d, weight)
            kk_bound = kk_maximum_weight_lower_bound(n, d, weight) if include_kk else 0
            dr_weight_bound = maximum_weight_dyachkov_rykov_lower_bound(
                n, d, weight
            )
            ruszinko_weight_bound = ruszinko_maximum_weight_lower_bound(
                n, d, weight
            )
            direct_weight_bound = maximum_weight_lower_bound(n, d, weight)
            iko_bound = inan_kairouz_ozgur_lower_bound(n, d, weight)
            sg_weight_bound = shangguan_ge_weight_aware_lower_bound(n, d, weight)
            combined_weight_bound = max(
                infinity, kk_bound, recursive_bound, direct_weight_bound,
                iko_bound, sg_weight_bound,
            )
            weight_bound = combined_weight_bound
            average_bound = combined_average_weight_lower_bound(
                n, d, weight, general
            )
            average_profile_bound = max(
                general,
                average_weight_profile_relaxation_lower_bound(n, d, weight),
            )
            best = max(
                furedi, ruszinko, dr, ks, bdr, kk_bound, recursive_bound,
                ruszinko_weight_bound, dr_weight_bound, direct_weight_bound,
                sg, iko_bound, sg_weight_bound,
            )
            best_average = average_bound
            prefix = {
                "EGH-general": "egh_general",
                "EGH-backtrack": "egh_backtrack",
                "Macula": "macula",
                "KS-B": "ks_b",
                "STW-SHF": "stw_shf",
                "KS-C": "ks_c",
                "VPR": "vpr",
                "GRV": "grv",
            }[construction]
            for name, value in (
                (f"{prefix}_m", m), (f"{prefix}_w", weight),
                (f"{prefix}_lower_recursive", recursive_bound),
                (f"{prefix}_lower_ruszinko_weight", ruszinko_weight_bound),
                (f"{prefix}_lower_dr_weight", dr_weight_bound),
                (f"{prefix}_lower_direct_weight", direct_weight_bound),
                (f"{prefix}_lower_iko", iko_bound),
                (f"{prefix}_lower_sg_weight", sg_weight_bound),
                (f"{prefix}_lower_combined_weight", combined_weight_bound),
                (f"{prefix}_lower_maximum", max(general, weight_bound)),
                (f"{prefix}_lower_weight", max(weight_bound, average_bound)),
                (f"{prefix}_lower_average", average_bound),
                (f"{prefix}_lower_average_profile", average_profile_bound),
                (f"{prefix}_gap", f"{m / best:.6f}"),
                (f"{prefix}_gap_average", f"{m / best_average:.6f}"),
                (f"{prefix}_summary", f"{m}/{weight}/{max(weight_bound, average_bound)}/{m / best:.3f}"),
            ):
                _text(row, name, value)
            if include_kk:
                _text(row, f"{prefix}_lower_kk", kk_bound)
                _text(row, f"{prefix}_kk_exact", int(kk_was_exact(d, weight)))
        inr_m = inr_parameters(n, d)
        _text(row, "inr_m", inr_m)
        _text(row, "inr_gap", f"{inr_m / general:.6f}")
        _text(row, "inr_summary", f"{inr_m}/-/{general}/{inr_m / general:.3f}")
        _text(row, "lower_furedi", furedi)
        _text(row, "lower_ruszinko", ruszinko)
        _text(row, "lower_dr", dr)
        _text(row, "lower_ks", ks)
        _text(row, "lower_sg", sg)
        _text(row, "lower_bdr", bdr)
        _text(row, "lower_infinity", infinity)
    return root


def write_xml(root: ET.Element, output: Path) -> None:
    ET.indent(root, space="  ")
    output.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--d", type=int, action="append", choices=DEFECTIVES,
        help="generate only this defective count (repeatable; default: all)",
    )
    parser.add_argument(
        "--include-kk", action="store_true",
        help="include Kruskal--Katona bounds in the XML and best-bound gaps",
    )
    args = parser.parse_args()
    for d in (args.d or DEFECTIVES):
        write_xml(
            measurement_xml(d, include_kk=args.include_kk),
            args.output_dir / f"cgt_gap_d{d}.xml",
        )


if __name__ == "__main__":
    main()
