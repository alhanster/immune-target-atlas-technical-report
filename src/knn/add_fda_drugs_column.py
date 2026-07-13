"""
Add an `FDA_drugs_targeting_nearest_gene` column to the KNN immune-target
candidate tables (both the Stim8hr and Stim48hr conditions).

For each candidate row we look up its `nearest_fda_target` gene in the
immune-system drug list and collect the FDA-approved (maxClinicalStage ==
"APPROVAL") drug names that target that gene. Drug names are kept verbatim
(salt-form variants included), de-duplicated, alphabetically sorted, and
joined with "; ".

For each input file, a new CSV is written next to it, leaving the original
untouched.
"""

import os
import warnings

# Silence harmless NumPy/pandas C-extension import warnings seen in this env.
warnings.filterwarnings("ignore")

import pandas as pd

# Paths anchored to this script's location so it runs from anywhere.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))   # src/knn
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))   # repo root
DERIVED = os.path.join(REPO_ROOT, "data", "derived")

# KNN candidate tables to process (one per culture condition).
KNN_CSVS = [
    os.path.join(DERIVED, "knn_immune_target_candidates_Stim8hr.csv"),
    os.path.join(DERIVED, "knn_immune_target_candidates_Stim48hr.csv"),
]
# Small committed subset (target -> approved drug names); produced by
# scripts/subset_approved_drugs.py from the raw Open Targets pull (see DATA.md).
DRUGS_CSV = os.path.join(
    REPO_ROOT, "data", "reference", "approved_target_drugs.csv"
)

NEW_COL = "FDA_drugs_targeting_nearest_gene"
TARGET_COL = "nearest_fda_target"


def build_gene_to_drugs(drugs_path):
    """Map each target gene -> "; "-joined string of FDA-approved drug names."""
    drugs = pd.read_csv(drugs_path)
    approved = drugs[drugs["maxClinicalStage"] == "APPROVAL"]

    gene_to_drugs = {}
    for gene, group in approved.groupby("Target (Gene)"):
        names = sorted(set(group["Drug Name"].dropna()))
        gene_to_drugs[gene] = "; ".join(names)
    return gene_to_drugs


def add_column(knn_path, gene_to_drugs):
    """Add the FDA-drugs column to one KNN file and write a new CSV alongside it."""
    knn = pd.read_csv(knn_path)

    values = knn[TARGET_COL].map(lambda g: gene_to_drugs.get(g, "")).fillna("")

    insert_at = knn.columns.get_loc(TARGET_COL) + 1
    if NEW_COL in knn.columns:
        knn = knn.drop(columns=[NEW_COL])
    knn.insert(insert_at, NEW_COL, values)

    root, ext = os.path.splitext(knn_path)
    out_path = f"{root}_with_fda_drugs{ext}"
    knn.to_csv(out_path, index=False)

    n_filled = (values != "").sum()
    print(f"Wrote {out_path}")
    print(f"  Rows: {len(knn)}  |  rows with >=1 FDA drug: {n_filled}")


def main():
    gene_to_drugs = build_gene_to_drugs(DRUGS_CSV)
    for knn_path in KNN_CSVS:
        add_column(knn_path, gene_to_drugs)


if __name__ == "__main__":
    main()
