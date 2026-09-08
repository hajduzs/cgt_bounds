# Reproduction package

This directory reproduces the numeric tables reported in the paper — every
lower bound, construction length, and gap-to-optimal ratio, computed from
first principles rather than transcribed. It is a deliberately small slice
of a much larger research codebase, containing only what the tables
actually depend on.

## Quick start

```
pip install .
./reproduce.sh
```

This writes the raw computed values to `output/measurements/*.xml`, the
LaTeX table fragments to `output/tables/*.tex`, and the LaTeX figure
fragment to `output/figures/*.tex`. Each is a drop-in `\input` target for
the paper:

| Generated file | Table / figure |
|---|---|
| `cgt_gap_d2_d3.tex`, `cgt_gap_d5_d10.tex` | Per-construction, per-`n` lower-bound and gap tables, for `d ∈ {2, 3, 5, 10}` |
| `cwc_small.tex` | Near-optimal small constant-weight codes (Brouwer-table comparison) |
| `construction_average_gaps.tex` | Mean gap-to-optimal per construction, per `d` |
| `ksb_weight_aware_bounds.tex` | Weight-aware vs. weight-independent lower bounds against the KS-B construction |

Reproduction is exact: running this script reproduces every table and the
figure byte-for-byte, verified against the full research repository this
package was extracted from.

## What's here, and why

Every lower bound and construction implemented here corresponds to a
specific theorem cited in the paper:

- **`src/cgt/bounds/`** — the lower-bound families compared against
  (Füredi, Ruszinkó, D'yachkov–Rykov, Kautz–Singleton, Shangguan–Ge,
  Inan–Kairouz–Özgür, the Kruskal–Katona refinement, and our own
  balanced-partition bound).
- **`src/cgt/constructions/`** — the explicit constructions compared
  (Eppstein–Goodrich–Hirschberg, Gargano–Rescigno–Vaccaro, Macula,
  Stinson–Trung–Wei).
- **`src/cgt/tables/`** — closed-form parameter formulas for the
  constructions that don't need explicit codewords (Kautz–Singleton's
  Reed–Solomon families, Porat–Rothschild, Indyk–Ngo–Rudra).
- **`src/cwc/bounds/brouwer.py`** — loads the tabulated small-length
  constant-weight-code bounds this comparison is checked against (see
  `data/`, below).
- **`src/core/`** — shared combinatorial primitives (binomial
  coefficients, the Johnson bound).
- **`scripts/`** — the five programs that turn the above into the
  measurement data, LaTeX tables, and figure; `xml_process.py` renders two
  of the tables from that data.

Nothing here is a registry, a CLI, or a search/optimization tool — the
source repository has all three, none of them are on the path from
"theorem" to "table," and stripping them out is most of what makes this
package small.

## Data

- **`data/code_lengths.json`** — previously-constructed constant-weight
  codes, used to report exact lengths for the recursive Kautz–Singleton
  family (KS-C) rather than only its closed-form bound.
- **`data/bounds/brouwer_*.txt`** — a parsed snapshot of Andries Brouwer's
  tables of best-known constant-weight-code bounds
  (<https://www.win.tue.nl/~aeb/codes/Andw.html>), used as an independent
  upper-bound reference in the small-code comparison. A dated, unmodified
  copy of the source pages this snapshot was parsed from — for provenance,
  not for re-parsing — can be provided alongside this package on request.

## Tests

```
pip install '.[test]'
pytest
```

The included tests check each bound and construction against its cited
theorem in isolation (hand-derived cases, known extremal values, and
cross-checks between independent formulas for the same quantity) — they
are the reviewable evidence that this code computes what the paper claims
it computes, independent of whether the generated tables happen to look
reasonable.
