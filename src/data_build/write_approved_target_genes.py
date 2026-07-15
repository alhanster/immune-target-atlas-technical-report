"""Derive the committed label file ``data/reference/approved_target_genes.txt``
from the gene rollup produced by ``build_approved_immune_ot.py``.

Input  : approved_immune_drugs_by_gene_ot.csv  (its ``gene_target`` column)
Output : data/reference/approved_target_genes.txt  (one gene symbol per line)

The output is the PU-model positive/anchor set consumed by src/model, src/knn,
and src/iei_enrichment (see DATA.md). Symbols are de-duplicated and sorted so the
committed file is stable across runs.

Usage:
    python src/data_build/write_approved_target_genes.py \
        --in  data/raw/approved_immune_drugs_by_gene_ot.csv \
        --out data/reference/approved_target_genes.txt
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_IN = REPO / "data" / "raw" / "approved_immune_drugs_by_gene_ot.csv"
DEFAULT_OUT = REPO / "data" / "reference" / "approved_target_genes.txt"
GENE_COL = "gene_target"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", type=Path, default=DEFAULT_IN)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    if not args.inp.exists():
        raise SystemExit(
            f"gene rollup not found: {args.inp}\n"
            "Regenerate it via src/data_build/build_approved_immune_ot.py (see DATA.md)."
        )

    genes: set[str] = set()
    with args.inp.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if GENE_COL not in (reader.fieldnames or []):
            raise SystemExit(
                f"column {GENE_COL!r} not in {args.inp.name}; "
                f"found {reader.fieldnames}"
            )
        for row in reader:
            # by_gene rows carry a single symbol; split defensively in case a
            # ';'-joined value is ever passed in.
            for sym in (row[GENE_COL] or "").split(";"):
                sym = sym.strip()
                if sym:
                    genes.add(sym)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as fh:
        for sym in sorted(genes):
            fh.write(sym + "\n")

    print(f"wrote {args.out} ({len(genes)} gene symbols) from {args.inp.name}")


if __name__ == "__main__":
    main()
