"""Subset the full gnomAD v4.1.1 constraint table to the few columns the
pipeline actually uses, so the repo can ship a small committed file instead of
the 278 MB original.

The full table (~110 columns) is only consumed by
``src/data_build/creating_full_gene_list.py``, which reads exactly:
``gene, gene_id, chromosome, lof.oe_ci.upper, mis.z_score``.
This script keeps only those columns and writes the committed subset.

Full source (see DATA.md):
    gnomAD v4.1.1 constraint metrics (canonical), downloadable from the gnomAD
    downloads page; place it at data/raw/ and point --input at it.

Usage:
    python scripts/subset_gnomad.py \
        --input data/raw/gnomad.v4.1.1.constraint_metrics.canonical.tsv \
        --output data/reference/gnomad_constraint_subset.tsv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
KEEP = ["gene", "gene_id", "chromosome", "lof.oe_ci.upper", "mis.z_score"]

DEFAULT_INPUT = REPO / "data" / "raw" / "gnomad.v4.1.1.constraint_metrics.canonical.tsv"
DEFAULT_OUTPUT = REPO / "data" / "reference" / "gnomad_constraint_subset.tsv"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    if not args.input.exists():
        raise SystemExit(
            f"gnomAD constraint table not found: {args.input}\n"
            "Download it (see DATA.md) into data/raw/ and pass --input."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(args.input, newline="") as fh_in, open(args.output, "w", newline="") as fh_out:
        reader = csv.DictReader(fh_in, delimiter="\t")
        missing = [c for c in KEEP if c not in reader.fieldnames]
        if missing:
            raise SystemExit(f"input is missing expected columns: {missing}")
        writer = csv.DictWriter(fh_out, fieldnames=KEEP, delimiter="\t")
        writer.writeheader()
        for row in reader:
            writer.writerow({c: row[c] for c in KEEP})
            n += 1

    size_mb = args.output.stat().st_size / 1e6
    print(f"wrote {args.output} ({n} rows, {size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
