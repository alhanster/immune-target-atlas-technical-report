"""PU-appropriate evaluation.

No plain accuracy / AUC-on-unlabeled. The primary metric is enrichment /
recall@k of HELD-OUT positives under group-aware CV: the bagging-PU is nested
inside each grouped fold and only ever sees the eval fold's genes at scoring
time, so a positive's family is never in both train and eval.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

from .config import Config
from .cv import make_folds
from .model import bagging_pu_predict_external


def cross_val_oof_scores(
    X: pd.DataFrame, y: np.ndarray, family: pd.Series, cfg: Config,
    model: str = "lgbm",
) -> np.ndarray:
    """Grouped-CV out-of-fold scores: every gene scored by models trained on
    other families. Bagging-PU is run inside each fold on that fold's train
    positives + a random equal draw of that fold's unlabeled genes.
    """
    oof = np.full(len(X), np.nan)
    folds = make_folds(family, y, cfg.n_folds, cfg.seed)
    for i, (tr, ev) in enumerate(folds):
        preds = bagging_pu_predict_external(
            X.iloc[tr], y[tr], X.iloc[ev], cfg, seed=cfg.seed + 1000 * (i + 1),
            model=model,
        )
        oof[ev] = preds
    assert not np.isnan(oof).any(), "some gene never appeared in an eval fold"
    return oof


def pu_metrics(y: np.ndarray, scores: np.ndarray, top_k_list) -> dict:
    """Recall@k, enrichment@k and PU average-precision on a ranked list.

    enrichment@k = (positive rate in top-k) / (overall positive rate); the
    fold-enrichment of held-out positives over a random ranking.
    """
    n = len(y)
    n_pos = int(y.sum())
    order = np.argsort(-scores, kind="mergesort")  # stable, high score first
    y_sorted = y[order]
    base_rate = n_pos / n
    out = {"n": n, "n_pos": n_pos, "base_rate": base_rate}
    for k in top_k_list:
        k = min(k, n)
        hits = int(y_sorted[:k].sum())
        recall = hits / n_pos if n_pos else 0.0
        precision = hits / k
        enrichment = (precision / base_rate) if base_rate else 0.0
        out[f"recall@{k}"] = recall
        out[f"hits@{k}"] = hits
        out[f"enrichment@{k}"] = enrichment
    # PU average precision (positives trusted; unlabeled treated as negatives ->
    # this is a conservative lower bound on true AP).
    out["pu_average_precision"] = float(average_precision_score(y, scores))
    return out
