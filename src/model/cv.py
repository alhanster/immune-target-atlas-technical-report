"""Group-aware cross-validation splitting.

Folds are grouped by gene family so that no family (hence no paralog cluster)
appears in both train and eval. Without this the model relearns "is-a-kinase /
GPCR / ion-channel" and CV enrichment is meaningless. Stratified on the label so
each fold carries a comparable share of the (rare) positives.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold


def make_folds(
    family: pd.Series, y: np.ndarray, n_folds: int, seed: int
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Return list of (train_idx, eval_idx) positional index arrays.

    family : per-gene family label (len == n_genes), aligned to y.
    y      : 0/1 label array.
    """
    groups = family.to_numpy()
    sgkf = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    folds = list(sgkf.split(np.zeros(len(y)), y, groups))
    # sanity: no family shared across the split
    for tr, ev in folds:
        assert not (set(groups[tr]) & set(groups[ev])), "family leaked across fold"
    return folds
