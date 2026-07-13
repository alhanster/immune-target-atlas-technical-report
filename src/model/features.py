"""Feature engineering and the per-column NA / imputation policy.

All 18,692 rows are kept. Missingness is handled per the task policy:
  - gwas_score          : already 0-imputed; add gwas_present indicator.
  - perturb-seq (9 cols): clean all-or-nothing block. impute0_indicator mode ->
                          fill 0 + ONE perturbseq_measured indicator; native mode
                          -> leave NaN for LightGBM, no indicator.
  - polarization_score  : impute 0 (or median) + polarization_measured indicator.
  - lof.oe_ci.upper / mis.z_score : genuinely MAR -> leave NaN (native), no indicator.

Derived features fuse effect+significance and capture dynamic response.
Family/druggability-proxy features are kept in a SEPARATE, ablatable group
(off by default) to avoid circularity.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .config import Config


def build_features(
    df: pd.DataFrame, cfg: Config, family: pd.Series | None = None
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """Build the model feature matrix.

    Returns (X, feature_groups, presence_flags) where feature_groups maps a
    group name -> list of columns (for ablation), and presence_flags holds the
    per-block measured indicators for the nominations output.
    """
    n = len(df)
    X = pd.DataFrame(index=df.index)
    groups: dict[str, list[str]] = {}

    # --- Genetics: gwas_score (ablatable via cfg.use_gwas) ----------------
    presence = pd.DataFrame(index=df.index)
    presence["gwas_present"] = (df["gwas_score"] != 0).astype(int)
    if cfg.use_gwas:
        X["gwas_score"] = df["gwas_score"].astype(float)
        X["gwas_present"] = presence["gwas_present"]
        groups["genetics"] = ["gwas_score", "gwas_present"]
    else:
        groups["genetics"] = []

    # --- IEI flag (inborn errors of immunity) -----------------------------
    X["IEI"] = df["IEI"].astype(int)
    groups.setdefault("prior_immune", []).append("IEI")

    # --- Constraint: MAR, leave NaN native, no indicator ------------------
    for c in config.CONSTRAINT_COLS:
        X[c] = pd.to_numeric(df[c], errors="coerce")
    groups["constraint"] = list(config.CONSTRAINT_COLS)

    # --- Perturb-seq block -------------------------------------------------
    pcols = config.PERTURBSEQ_COLS
    measured = df[pcols].notna().any(axis=1).astype(int)
    presence["perturbseq_measured"] = measured
    keep_pcols = list(pcols)
    if cfg.drop_coef_se:
        keep_pcols = [c for c in pcols if c not in config.COEF_SE_COLS]

    if cfg.assay_na_mode == "impute0_indicator":
        for c in keep_pcols:
            X[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
        X["perturbseq_measured"] = measured
        groups["perturbseq"] = keep_pcols + ["perturbseq_measured"]
    else:  # native: leave NaN, no indicator
        for c in keep_pcols:
            X[c] = pd.to_numeric(df[c], errors="coerce")
        groups["perturbseq"] = list(keep_pcols)

    # --- Derived perturb-seq features (computed from raw, NaN-aware) -------
    slp = df[config.SIGNED_LOG10P_COLS].apply(pd.to_numeric, errors="coerce")
    max_abs_slp = slp.abs().max(axis=1)
    dyn_slp = pd.to_numeric(df["stim48hr_signed_log10p"], errors="coerce") - \
        pd.to_numeric(df["rest_signed_log10p"], errors="coerce")
    dyn_coef = pd.to_numeric(df["stim48hr_coef_beta"], errors="coerce") - \
        pd.to_numeric(df["rest_coef_beta"], errors="coerce")
    responsive_any = (slp.abs() >= 1.0).any(axis=1).astype(float)  # |signed_log10p|>=1 ~ p<=0.1

    if cfg.assay_na_mode == "impute0_indicator":
        max_abs_slp = max_abs_slp.fillna(0.0)
        dyn_slp = dyn_slp.fillna(0.0)
        dyn_coef = dyn_coef.fillna(0.0)
        responsive_any = responsive_any.fillna(0.0)
    X["max_abs_signed_log10p"] = max_abs_slp
    X["dyn_signed_log10p_48_rest"] = dyn_slp
    X["dyn_coef_48_rest"] = dyn_coef
    X["responsive_any"] = responsive_any
    groups["perturbseq_derived"] = [
        "max_abs_signed_log10p", "dyn_signed_log10p_48_rest",
        "dyn_coef_48_rest", "responsive_any",
    ]

    # --- Polarization: 80% NA -> impute (0 or median) + indicator ---------
    pol = pd.to_numeric(df["polarization_score"], errors="coerce")
    presence["polarization_measured"] = pol.notna().astype(int)
    fill = 0.0 if cfg.polarization_impute == "zero" else float(pol.median())
    X["polarization_score"] = pol.fillna(fill)
    X["polarization_measured"] = presence["polarization_measured"]
    groups["polarization"] = ["polarization_score", "polarization_measured"]

    # --- Family / druggability-proxy group (ablatable, OFF by default) ----
    if cfg.use_druggability and family is not None:
        fam = family.reindex(df["gene"].values).reset_index(drop=True)
        fam.index = df.index
        fam_size = fam.map(fam.value_counts()).astype(float)
        X["gene_family_size"] = fam_size
        X["in_hgnc_group"] = fam.str.startswith("HGNC:").astype(int)
        groups["family_druggability"] = ["gene_family_size", "in_hgnc_group"]

    assert len(X) == n, "row count changed — rows must never be dropped"
    return X, groups, presence


def active_feature_columns(groups: dict, use_druggability: bool) -> list[str]:
    """Flatten feature groups into the ordered column list actually used."""
    cols: list[str] = []
    for name, gcols in groups.items():
        if name == "family_druggability" and not use_druggability:
            continue
        cols.extend(gcols)
    # de-dup preserving order
    seen, out = set(), []
    for c in cols:
        if c not in seen:
            seen.add(c); out.append(c)
    return out
