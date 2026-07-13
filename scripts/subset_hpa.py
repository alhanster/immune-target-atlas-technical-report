"""Reduce the Human Protein Atlas immune-cell RNA table to per-gene max nTPM,
so the repo ships a tiny committed file instead of the 17 MB long-format original.

``src/iei_enrichment/make_enrichment_table.py`` only ever uses, per gene, the
maximum nTPM across immune cell types (to define the "immune-expressed" universe
at nTPM >= 1 / >= 5). This collapses the ~383k long rows to one row per gene while
keeping the exact columns the enrichment script reads (``Gene``, ``Gene name``,
``nTPM``), so its ``groupby(...).max()`` stays correct (idempotent) on the subset.

Full source (see DATA.md): Human Protein Atlas "RNA immune cell" TSV
(rna_immune_cell.tsv), downloadable from proteinatlas.org; place under data/raw/.

Usage:
    python scripts/subset_hpa.py \
        --input data/raw/rna_immune_cell.tsv \
        --output data/reference/hpa_immune_max_ntpm.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO / "data" / "raw" / "rna_immune_cell.tsv"
DEFAULT_OUTPUT = REPO / "data" / "reference" / "hpa_immune_max_ntpm.tsv"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    if not args.input.exists():
        raise SystemExit(
            f"HPA table not found: {args.input}\n"
            "Download rna_immune_cell.tsv (see DATA.md) into data/raw/ and pass --input."
        )

    df = pd.read_csv(args.input, sep="\t")
    reduced = (
        df.groupby(["Gene", "Gene name"])["nTPM"].max().reset_index()
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    reduced.to_csv(args.output, sep="\t", index=False)
    size_kb = args.output.stat().st_size / 1e3
    print(f"wrote {args.output} ({len(reduced)} genes, {size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
