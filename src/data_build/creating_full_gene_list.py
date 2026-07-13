#!/usr/bin/env python3
"""Filter the gnomAD v4.1.1 canonical constraint metrics table.

Keeps rows where:
  1. gene_id is an Ensembl gene ID (starts with "ENSG")
  2. the gene is NOT on a sex chromosome (chrX / chrY)
  3. the gene name is present (not "NA")

Then drops every row whose gene *symbol* is duplicated: a handful of
symbols map to two distinct Ensembl gene IDs, and with no rule to pick a
"correct" transcript we remove all copies so each symbol is unique.

Finally, harmonizes gene symbols to current HGNC-approved symbols via the
shared `gene_name_utils` module: outdated "previous" symbols are renamed,
and a deprecated symbol whose current name already exists in the list is
dropped so the correct existing row is kept (no duplicate). The rename map
is derived dynamically from the HGNC complete set.

Also adds an `IEI` column: 1 if the gene is in the curated IUIS
inborn-errors-of-immunity panel (IEI_gene_list.csv), else 0.

Finally, left-joins three external per-gene score tables by gene name
(each harmonized to the same HGNC release first): the GWAS gene-level
score (gwas_gene_scores_*.csv), the regulator-burden wide scores
(regulator_burden_scores_wide.csv), and the polarization score
(polarization_gene_scores.csv). Genes missing from the burden/polarization
assays get NA (not measured); genes missing a GWAS score get 0 (no
credible-set association signal, i.e. a real zero).

Usage:
    python creating_full_gene_list.py \
        [input.tsv] [output.tsv]

Defaults to the canonical file in this directory and writes
`gnomad.v4.1.1.constraint_metrics.canonical.filtered.tsv` alongside it.
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent        # src/data_build
REPO = HERE.parents[1]                         # repo root
sys.path.insert(0, str(REPO / "src" / "shared"))   # so `import gene_name_utils` resolves
from gene_name_utils import load_hgnc, harmonize

# gnomAD subset (committed; produced by scripts/subset_gnomad.py — see DATA.md).
DEFAULT_INPUT = REPO / "data" / "reference" / "gnomad_constraint_subset.tsv"
DEFAULT_OUTPUT = REPO / "data" / "derived" / "full_gene_list.tsv"
DEFAULT_IEI = REPO / "data" / "reference" / "IEI_gene_list.csv"

# Sex chromosomes to drop (the table uses the "chr" prefix, e.g. "chrX").
SEX_CHROMOSOMES = {"chrX", "chrY", "X", "Y"}

# Columns kept from the constraint table, in this order (join columns appended).
OUTPUT_COLUMNS = ["gene", "gene_id", "lof.oe_ci.upper", "mis.z_score", "IEI"]

# External per-gene score tables left-joined onto the list by gene name, in order.
# Each is harmonized to current HGNC symbols before joining. Tuple fields:
#   (path, value_cols, fill)
#   value_cols=None -> "all columns except gene".
#   fill=None       -> unmatched genes stay NA (score not measured for that gene).
#   fill=0          -> unmatched genes get 0. Used for the GWAS score: absence of
#                      a credible-set association means no signal, i.e. a real 0
#                      (not "unmeasured"), unlike the burden/polarization assays.
DERIVED = REPO / "data" / "derived"
JOIN_SOURCES = [
    (DERIVED / "gwas_gene_scores.csv",               ["gwas_score"],         0),
    (DERIVED / "regulator_burden_scores_wide.csv",   None,                   None),
    (DERIVED / "polarization_gene_scores.csv",       ["polarization_score"], None),
]


def _load_scores(path: Path, hgnc: dict, value_cols):
    """Load a per-gene score CSV, harmonize its gene symbols to current HGNC,
    and return ['gene', *value_cols] ready for a left join (value_cols=None
    keeps every column except 'gene')."""
    src = pd.read_csv(path)
    src, _ = harmonize(src, hgnc, gene_col="gene", verbose=False)
    if value_cols is None:
        value_cols = [c for c in src.columns if c != "gene"]
    return src[["gene", *value_cols]]


def filter_constraint_table(
    input_path: Path, output_path: Path, iei_path: Path
) -> None:
    # keep_default_na=False so the literal string "NA" is preserved as text
    # (we filter it explicitly) rather than being coerced to a float NaN.
    df = pd.read_csv(input_path, sep="\t", dtype=str, keep_default_na=False)
    n_start = len(df)

    # 1. gene_id must be an Ensembl gene ID.
    mask_ensg = df["gene_id"].str.startswith("ENSG")

    # 2. drop sex chromosomes.
    mask_autosomal = ~df["chromosome"].isin(SEX_CHROMOSOMES)

    # 3. drop missing gene names.
    mask_named = df["gene"].str.strip().ne("") & df["gene"].ne("NA")

    filtered = df[mask_ensg & mask_autosomal & mask_named]
    n_filtered = len(filtered)

    # Drop every row whose gene symbol appears more than once (keep=False
    # marks all copies, so both members of each collision are removed).
    dup_mask = filtered["gene"].duplicated(keep=False)
    dup_names = sorted(filtered.loc[dup_mask, "gene"].unique())
    deduped = filtered[~dup_mask]

    # Harmonize to current HGNC-approved symbols: rename outdated "previous"
    # symbols, and drop a deprecated symbol whose current name already exists
    # (keeping the correct existing row). Dynamic map from the HGNC complete set.
    hgnc = load_hgnc()                                    # loaded once, reused below
    deduped, hginfo = harmonize(deduped, hgnc, gene_col="gene", verbose=False)
    n_renamed = hginfo["renamed"]
    hgnc_dropped = hginfo["dropped"]

    # Flag membership in the curated IUIS inborn-errors-of-immunity panel.
    # Both sides use current HGNC-approved symbols (symbols were harmonized just
    # above), so a direct symbol match is correct.
    iei_symbols = set(
        pd.read_csv(iei_path, dtype=str, keep_default_na=False)["Gene"]
    )
    deduped = deduped.assign(IEI=deduped["gene"].isin(iei_symbols).astype(int))

    # Keep the constraint columns, then left-join the external score tables by
    # gene name (harmonized to the same HGNC release). Genes missing from a
    # source stay NA, except sources with a fill value (GWAS -> 0).
    selected = deduped[OUTPUT_COLUMNS]
    join_counts = {}
    for path, value_cols, fill in JOIN_SOURCES:
        scores = _load_scores(path, hgnc, value_cols)
        joined_cols = [c for c in scores.columns if c != "gene"]
        selected = selected.merge(scores, on="gene", how="left")
        # Count matches on the first joined column before any fill is applied.
        join_counts[path.name] = int(selected[joined_cols[0]].notna().sum())
        if fill is not None:
            selected[joined_cols] = selected[joined_cols].fillna(fill)

    # na_rep so unmatched (NaN) score cells render as the literal "NA" (matching
    # the constraint table's existing string NAs).
    selected.to_csv(output_path, sep="\t", index=False, na_rep="NA")

    print(f"Input rows:              {n_start:,}")
    print(f"  gene_id starts ENSG:   {int(mask_ensg.sum()):,}")
    print(f"  not on sex chromosome: {int(mask_autosomal.sum()):,}")
    print(f"  gene name present:     {int(mask_named.sum()):,}")
    print(f"Rows kept (all filters): {n_filtered:,}")
    print(f"Duplicated symbols dropped ({len(dup_names)}): {', '.join(dup_names)}")
    print(f"  rows removed:          {int(dup_mask.sum()):,}")
    print(f"Rows kept (unique gene): {n_filtered - int(dup_mask.sum()):,}")
    print(f"HGNC deprecated dropped: {len(hgnc_dropped)} ({', '.join(hgnc_dropped)})")
    print(f"HGNC symbols renamed:    {n_renamed:,}")
    print(f"Rows kept (final):       {len(deduped):,}")
    print(f"IEI genes flagged (=1):  {int(selected['IEI'].sum())} / {len(iei_symbols)} in panel")
    for name, n in join_counts.items():
        print(f"Joined {name}: {n:,} genes matched (rest NA)")
    print(f"Output shape:            {selected.shape[0]:,} rows x {selected.shape[1]} cols")
    print(f"Wrote: {output_path}")


def main() -> None:
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INPUT
    output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUTPUT

    if not input_path.exists():
        sys.exit(f"Input file not found: {input_path}")
    if not DEFAULT_IEI.exists():
        sys.exit(f"IEI gene list not found: {DEFAULT_IEI}")
    for path, *_ in JOIN_SOURCES:
        if not path.exists():
            sys.exit(f"Join source not found: {path}")

    filter_constraint_table(input_path, output_path, DEFAULT_IEI)


if __name__ == "__main__":
    main()
