#!/usr/bin/env python3
"""Reshape the genome-wide regulator-burden scores to one row per gene.

Reads the long-format `regulator_burden_scores_all.csv` (one row per
gene-condition) from Part 2 and pivots it wide: for each of the three culture
conditions (Rest, Stim8hr, Stim48hr) it emits `coef_beta`, `se_beta`, and
`signed_log10p` columns, so each gene occupies a single row with no duplicates.

Also harmonizes gene symbols to current HGNC-approved symbols and runs a
gene-name formatting check, both via the shared `gene_name_utils` module (which
derives the rename map dynamically from the HGNC complete set).

Output: `regulator_burden_scores_wide.csv` in this directory.
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent        # src/data_build
REPO = HERE.parents[1]                         # repo root
sys.path.insert(0, str(REPO / "src" / "shared"))   # so `import gene_name_utils` resolves
from gene_name_utils import load_hgnc, harmonize, check_format
INPUT = REPO / "data" / "derived" / "regulator_burden_scores_all.csv"
OUTPUT = REPO / "data" / "derived" / "regulator_burden_scores_wide.csv"

# Metrics carried across (source name -> output suffix).
METRICS = {"coef_beta": "coef_beta", "se_beta": "se_beta", "signed_log10p": "signed_log10p"}
# Conditions in output order, with their lowercased column prefix.
CONDITIONS = [("Rest", "rest"), ("Stim8hr", "stim8hr"), ("Stim48hr", "stim48hr")]



def pivot_wide(df: pd.DataFrame) -> pd.DataFrame:
    """One row per gene; a {cond}_{metric} column for every condition/metric."""
    sub = df[["gene", "condition", *METRICS]]
    wide = sub.pivot(index="gene", columns="condition", values=list(METRICS))

    # Flatten ("coef_beta", "Rest") -> "rest_coef_beta", ordered by condition.
    cols = []
    for cond, prefix in CONDITIONS:
        for src, suffix in METRICS.items():
            cols.append((f"{prefix}_{suffix}", (src, cond)))
    out = pd.DataFrame({name: wide[key] for name, key in cols})
    return out.reset_index()


def main() -> None:
    if not INPUT.exists():
        raise SystemExit(f"Input not found: {INPUT}")
    df = pd.read_csv(INPUT)

    wide = pivot_wide(df)
    wide, _ = harmonize(wide, load_hgnc(), gene_col="gene")
    check_format(wide["gene"])

    wide.to_csv(OUTPUT, index=False)
    print(f"\nwrote {OUTPUT}  shape={wide.shape}")
    print("columns:", list(wide.columns))


if __name__ == "__main__":
    main()
