#!/usr/bin/env python3
"""
Gene-level polarization score for CD4+ T-cell Perturb-seq.

polarization_score = (z_Stim8hr + z_Stim48hr) / sqrt(2)          # Stouffer combination
    where z_condition = coef_mean / coef_sem  on the Th1/Th2 ("ota") signature.

Scores every regulator measured in both stimulation timepoints (Rest optional,
kept only as an annotation) -- not restricted to any gene panel.

Input:  polarization_prediction_condition_comparison_regulator_coefficients.csv
        (metadata/suppl_tables/ in github.com/emdann/GWT_perturbseq_analysis_2025)
Output: polarization_score.csv

The diagnostic figure is drawn separately by make_polarization_score_figure.R,
which reads polarization_score.csv -- this script does not produce figures.

Usage:   python polarization_score.py [coefficients.csv]
         (defaults to ../inputs/polarization_prediction_condition_comparison_regulator_coefficients.csv)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# Paths anchored to this script's location so it runs from any cwd. The vendored
# coefficient table is vendored into this stage's inputs/ (see inputs/SOURCE.md).
HERE = Path(__file__).resolve().parent        # src/regulator_burden
REPO = HERE.parents[1]                         # repo root
DEFAULT_IN = (HERE / "inputs"
              / "polarization_prediction_condition_comparison_regulator_coefficients.csv")
IN_CSV  = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_IN
OUT_CSV = REPO / "data" / "derived" / "polarization_score.csv"

# Shared gene-name utilities live in src/shared/.
sys.path.insert(0, str(REPO / "src" / "shared"))
from gene_name_utils import load_hgnc, harmonize, check_format   # noqa: E402

# Slim, HGNC-harmonized deliverable (gene + polarization_score) for downstream use.
SLIM_OUT = REPO / "data" / "derived" / "polarization_gene_scores.csv"


# --------------------------------------------------------------------------- #
# 1. Load coefficients and build confidence-weighted z per condition
# --------------------------------------------------------------------------- #
def load_z(path):
    pc  = pd.read_csv(path, index_col=0)
    ota = pc[pc["signature"] == "ota"].copy()          # Th1/Th2 polarization signature
    conds = ["Rest", "Stim8hr", "Stim48hr"]
    coef = ota.pivot_table(index="regulator", columns="celltype",
                           values="coef_mean", aggfunc="first").reindex(columns=conds)
    sem  = ota.pivot_table(index="regulator", columns="celltype",
                           values="coef_sem",  aggfunc="first").reindex(columns=conds)
    z = coef / sem
    known = set(ota[ota["known_regulators"] == True]["regulator"])
    rtype = (ota.dropna(subset=["regulator_type"])
                .drop_duplicates("regulator")
                .set_index("regulator")["regulator_type"])
    return coef, z, known, rtype


# --------------------------------------------------------------------------- #
# 2. Compute the score
# --------------------------------------------------------------------------- #
def compute_score(coef, z, known, rtype):
    S = pd.DataFrame(index=coef.index)
    S["z_Stim8hr"]  = z["Stim8hr"]
    S["z_Stim48hr"] = z["Stim48hr"]
    S["z_Rest"]     = z["Rest"]
    # Require the two stimulation timepoints the score is built from (a valid
    # 2-term Stouffer combination); Rest is optional (annotation only, so its
    # z may be NaN for genes measured in stim but not rest).
    S = S[coef[["Stim8hr", "Stim48hr"]].notna().all(axis=1)].copy()

    # THE SCORE: Stouffer combination of the two stimulation timepoints
    S["polarization_score"] = (S["z_Stim8hr"] + S["z_Stim48hr"]) / np.sqrt(2)
    S["abs_score"]   = S["polarization_score"].abs()
    S["direction"]   = np.where(S["polarization_score"] > 0,
                                "toward_signature_pos", "toward_signature_neg")
    S["signs_agree"] = np.sign(S["z_Stim8hr"]) == np.sign(S["z_Stim48hr"])
    S["z_Rest_abs"]  = S["z_Rest"].abs()               # retained for window filtering

    S["coef_Stim8hr"]    = coef["Stim8hr"]
    S["coef_Stim48hr"]   = coef["Stim48hr"]
    S["known_regulator"] = S.index.isin(known)
    S["regulator_type"]  = rtype.reindex(S.index)

    S = S.sort_values("abs_score", ascending=False)
    S["rank"] = range(1, len(S) + 1)
    return S


# --------------------------------------------------------------------------- #
# 3. Calibration against study-annotated Th1/Th2 regulators
# --------------------------------------------------------------------------- #
def calibrate(S):
    kn = S[S["known_regulator"]]["abs_score"]
    bg = S[~S["known_regulator"]]["abs_score"]
    _, p = mannwhitneyu(kn, bg, alternative="greater")
    return p, len(kn)


# --------------------------------------------------------------------------- #
def main():
    coef, z, known, rtype = load_z(IN_CSV)
    S = compute_score(coef, z, known, rtype)
    p, n_known = calibrate(S)
    print(f"genes scored: {len(S)}")
    print(f"score range: [{S['polarization_score'].min():.1f}, {S['polarization_score'].max():.1f}] (signed z)")
    print(f"calibration: known Th1/Th2 regulators (n={n_known}) vs background, |score| MWU p = {p:.2e}")

    cols = ["gene", "rank", "polarization_score", "abs_score", "direction", "signs_agree",
            "z_Stim8hr", "z_Stim48hr", "z_Rest", "z_Rest_abs",
            "coef_Stim8hr", "coef_Stim48hr", "known_regulator", "regulator_type"]
    S.reset_index(names="gene")[cols].to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}")

    # Slim, HGNC-harmonized deliverable: gene + polarization_score only.
    slim = S.reset_index(names="gene")[["gene", "polarization_score"]]
    slim, _ = harmonize(slim, load_hgnc(), gene_col="gene")
    check_format(slim["gene"])
    slim.to_csv(SLIM_OUT, index=False)
    print(f"wrote {SLIM_OUT}")


if __name__ == "__main__":
    main()
