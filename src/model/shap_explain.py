"""Per-gene SHAP explanations for the top-N nominated unlabeled genes.

The bagging-PU score is the mean over base models, so the ensemble SHAP value is
the mean of the base models' SHAP values. We average TreeExplainer contributions
across the bags on just the top-N rows (cheap) and tag each gene by whether its
score is driven by genetics+Perturb-seq evidence or by the family/druggability
group.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import shap

GENETICS_PERTURB_PREFIXES = (
    "gwas", "IEI", "max_abs_signed_log10p", "dyn_", "responsive_any",
    "rest_", "stim8hr_", "stim48hr_", "perturbseq_measured",
)
FAMILY_DRUGGABILITY = ("gene_family_size", "in_hgnc_group")


def explain_top_genes(
    bags: list, X: pd.DataFrame, feat_cols: list[str], row_idx: np.ndarray,
) -> pd.DataFrame:
    """Return a tidy DataFrame of mean SHAP values (over bags) for the given rows.

    One row per (gene position, feature); plus a per-gene driver tag is computed
    by the caller from these values.
    """
    Xr = X.iloc[row_idx][feat_cols]
    shap_sum = np.zeros((len(row_idx), len(feat_cols)))
    for clf in bags:
        booster = getattr(clf, "booster_", None)
        expl = shap.TreeExplainer(booster if booster is not None else clf)
        vals = expl.shap_values(Xr)
        if isinstance(vals, list):  # older shap: [class0, class1]
            vals = vals[1]
        shap_sum += np.asarray(vals)
    shap_mean = shap_sum / len(bags)
    return pd.DataFrame(shap_mean, columns=feat_cols, index=row_idx)


def driver_tag(shap_row: pd.Series) -> str:
    """Classify a gene as genetics+perturbseq-driven vs family/druggability-driven
    using the share of positive SHAP mass."""
    pos = shap_row.clip(lower=0)
    total = pos.sum()
    if total <= 0:
        return "weak/mixed"
    gp = sum(v for f, v in pos.items()
             if f.startswith(GENETICS_PERTURB_PREFIXES))
    fam = sum(v for f, v in pos.items() if f in FAMILY_DRUGGABILITY)
    gp_frac, fam_frac = gp / total, fam / total
    if fam_frac > 0.5:
        return "family/druggability-driven"
    if gp_frac >= 0.5:
        return "genetics+perturbseq-driven"
    return "mixed"


def top_feature_string(shap_row: pd.Series, k: int = 3) -> str:
    """Human-readable top-k contributing features (signed)."""
    top = shap_row.reindex(shap_row.abs().sort_values(ascending=False).index).head(k)
    return "; ".join(f"{f}{'+' if v >= 0 else '-'}{abs(v):.3f}" for f, v in top.items())
