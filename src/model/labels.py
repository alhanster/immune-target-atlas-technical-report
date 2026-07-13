"""Label construction and the gwas_score provenance guardrail.

Two positive sets are always built (guardrail #2):
  a) 'all'    — every approved FDA target present in the gene list.
  b) 'immune' — approved targets that ALSO have gwas_score>0 OR IEI==1 (the
                primary/honest immune-restricted target).
"""
from __future__ import annotations

import pandas as pd

# ---------------------------------------------------------------------------
# Guardrail #1: gwas_score provenance.
#
# Verified from `src/gwas/gwas_gene_scores.py`: the score is the
# Open Targets `gwas_credible_sets` DATASOURCE score (the GWAS component of the
# `genetic_association` datatype), pulled via an associatedTargets query with
# `datasources: [{ id: "gwas_credible_sets", required: true }]` and read from
# `datasourceScores` where id == "gwas_credible_sets". It is NOT the overall
# association score, so it does NOT fold in the ChEMBL known-drug datasource and
# therefore does NOT leak approval labels.
# ---------------------------------------------------------------------------
GWAS_PROVENANCE_CONFIRMED = True
GWAS_PROVENANCE_NOTE = (
    "gwas_score = Open Targets 'gwas_credible_sets' datasource score for "
    "MONDO:0005046 (genetics/GWAS component of genetic_association), NOT the "
    "overall association score. Confirmed from src/gwas/gwas_gene_scores.py "
    "(datasources filter + datasourceScores read). No ChEMBL/known-drug leakage."
)


def gwas_provenance_warning() -> str:
    """Return the loud provenance banner to print at run start."""
    if GWAS_PROVENANCE_CONFIRMED:
        return (
            "[GUARDRAIL 1 — gwas_score provenance: CONFIRMED SAFE]\n"
            f"  {GWAS_PROVENANCE_NOTE}\n"
            "  Ablation with/without gwas_score is still run for transparency."
        )
    return (
        "[GUARDRAIL 1 — gwas_score provenance: *** UNCONFIRMED ***]\n"
        "  Could not confirm gwas_score is the genetics/GWAS-datasource score.\n"
        "  If it is the OVERALL association score it folds in ChEMBL known-drug\n"
        "  evidence and LEAKS the labels. Running BOTH with and without it."
    )


def build_labels(df: pd.DataFrame, fda: set[str]) -> pd.DataFrame:
    """Add label columns to df (returns the same frame).

    Adds:
      label_all      — 1 if gene is an approved FDA target.
      label_immune   — 1 if approved AND (gwas_score>0 OR IEI==1).
      gwas_present   — gwas_score != 0 indicator.
    """
    df = df.copy()
    df["label_all"] = df["gene"].isin(fda).astype(int)
    df["gwas_present"] = (df["gwas_score"] != 0).astype(int)
    immune_evidence = (df["gwas_score"] > 0) | (df["IEI"] == 1)
    df["label_immune"] = (df["label_all"].astype(bool) & immune_evidence).astype(int)
    return df


def label_summary(df: pd.DataFrame, fda: set[str]) -> dict:
    n_fda_file = len(fda)
    n_in_list = int(df["gene"].isin(fda).sum())
    return {
        "n_fda_file": n_fda_file,
        "n_fda_in_gene_list": n_in_list,
        "n_fda_missing_from_list": n_fda_file - n_in_list,
        "positives_all": int(df["label_all"].sum()),
        "positives_immune": int(df["label_immune"].sum()),
        "n_genes": len(df),
    }
