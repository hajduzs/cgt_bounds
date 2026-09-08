#!/usr/bin/env bash
# Regenerates every numeric table and figure reported in the paper from
# source.
#
# Usage:  ./reproduce.sh [output-dir]
# Output: output-dir/measurements/*.xml   (raw computed values)
#         output-dir/tables/*.tex         (the LaTeX table fragments)
#         output-dir/figures/*.tex        (the PGFPlots figure fragment)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$ROOT/output}"
MEASUREMENTS="$OUT/measurements"
TABLES="$OUT/tables"
FIGURES="$OUT/figures"

export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
PYTHON="${PYTHON:-python3}"

mkdir -p "$MEASUREMENTS" "$TABLES" "$FIGURES"

echo "[1/5] Computing lower bounds and construction lengths for d in {2,3,5,10}..."
for d in 2 3 5 10; do
    "$PYTHON" "$ROOT/scripts/generate_article_gap_xml.py" \
        --output-dir "$MEASUREMENTS" --d "$d"
done

echo "[2/5] Computing the small constant-weight-code comparison..."
"$PYTHON" "$ROOT/scripts/generate_article_cwc_xml.py" \
    --output "$MEASUREMENTS/cwc_small.xml"

echo "[3/5] Rendering the detailed per-construction gap tables..."
"$PYTHON" "$ROOT/scripts/generate_article_detailed_tables.py" \
    --xml-dir "$MEASUREMENTS" --output-dir "$TABLES"

echo "[4/5] Rendering the KS-B weight-aware bounds figure..."
"$PYTHON" "$ROOT/scripts/generate_ksb_weight_aware_bounds_figure.py" \
    --xml-dir "$MEASUREMENTS" --output "$FIGURES/ksb_weight_aware_bounds.tex"

echo "[5/5] Rendering the small-CWC and mean-gap summary tables..."
"$PYTHON" "$ROOT/xml_process.py" "$MEASUREMENTS/cwc_small.xml" \
    -aggregate measurement -x row_id -hidex \
    -columns avg:d avg:m avg:weight avg:distance avg:size avg:lower_weight avg:gap \
    -columnnames d m w D n '$\mathbf L_w$' 'GAP$_w$' \
    -latex-column-format 'rrrrr|rr' -latex-double-header-rule \
    -latex -precision 3 -outfile "$TABLES/cwc_small.tex" >/dev/null
"$PYTHON" "$ROOT/scripts/generate_article_average_table.py" \
    --xml-dir "$MEASUREMENTS" --output "$TABLES/construction_average_gaps.tex"

echo
echo "Done. Table fragments written to: $TABLES"
echo "  cgt_gap_d2_d3.tex, cgt_gap_d5_d10.tex   -- detailed per-n, per-construction gap tables"
echo "  cwc_small.tex                           -- near-optimal small constant-weight codes"
echo "  construction_average_gaps.tex           -- mean GAP_best per construction, per d"
echo "Figure fragment written to: $FIGURES/ksb_weight_aware_bounds.tex"
echo
echo "These \\input directly into the paper; see README.md for the exact mapping."
