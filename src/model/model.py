"""Bagging-PU model and the elastic-net logistic baseline.

Bagging-PU (Mordelet & Vert): train N base learners, each on ALL positives plus
a random equal-size draw of unlabeled genes treated as pseudo-negatives; average
the base scores. For the FINAL ranking, each gene is scored only on the models
where it was out-of-bag (never drawn as a pseudo-negative) — an OOB average.
Positives are in every bag, so they are scored in-bag and flagged as such.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

from .config import Config


def _lgbm(cfg: Config, seed: int) -> LGBMClassifier:
    return LGBMClassifier(
        n_estimators=cfg.n_estimators,
        num_leaves=cfg.num_leaves,
        min_child_samples=cfg.min_child_samples,
        learning_rate=cfg.learning_rate,
        subsample=cfg.subsample,
        subsample_freq=1,
        colsample_bytree=cfg.colsample_bytree,
        random_state=seed,
        n_jobs=-1,
        verbose=-1,
    )


def _elasticnet() -> Pipeline:
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5,
            C=1.0, max_iter=2000, n_jobs=-1)),
    ])


def _draw_unlabeled(rng, unlabeled_idx: np.ndarray, n_pos: int) -> np.ndarray:
    """Draw n_pos pseudo-negatives from the unlabeled pool (without replacement
    if the pool is large enough, else with replacement)."""
    replace = len(unlabeled_idx) < n_pos
    return rng.choice(unlabeled_idx, size=n_pos, replace=replace)


def bagging_pu_predict_external(
    X_train: pd.DataFrame, y_train: np.ndarray, X_eval: pd.DataFrame,
    cfg: Config, seed: int, model: str = "lgbm",
) -> np.ndarray:
    """Train bagging-PU on (X_train, y_train) and average predictions on X_eval.

    Used for grouped-CV: X_eval is a disjoint family fold, so every bag scores it
    (no OOB bookkeeping needed).
    """
    rng = np.random.default_rng(seed)
    pos_idx = np.where(y_train == 1)[0]
    unl_idx = np.where(y_train == 0)[0]
    n_pos = len(pos_idx)
    preds = np.zeros(len(X_eval))
    Xtr = X_train.to_numpy()
    Xev = X_eval.to_numpy()
    for b in range(cfg.n_bags):
        neg = _draw_unlabeled(rng, unl_idx, n_pos)
        rows = np.concatenate([pos_idx, neg])
        yb = np.concatenate([np.ones(n_pos), np.zeros(len(neg))])
        clf = _lgbm(cfg, seed + b) if model == "lgbm" else _elasticnet()
        clf.fit(Xtr[rows], yb)
        preds += clf.predict_proba(Xev)[:, 1]
    return preds / cfg.n_bags


def bagging_pu_full(
    X: pd.DataFrame, y: np.ndarray, cfg: Config, seed: int, model: str = "lgbm",
):
    """Train bagging-PU on all data; return OOB-averaged scores for every gene.

    Returns (scores, in_bag_count, oob_count, bags) where scores use OOB
    averaging for unlabeled genes and in-bag averaging for positives (always
    in-bag). `bags` is the list of fitted models (for SHAP).
    """
    rng = np.random.default_rng(seed)
    pos_idx = np.where(y == 1)[0]
    unl_idx = np.where(y == 0)[0]
    n_pos = len(pos_idx)
    n = len(X)
    Xn = X.to_numpy()

    oob_sum = np.zeros(n)
    oob_cnt = np.zeros(n)
    inbag_sum = np.zeros(n)
    inbag_cnt = np.zeros(n)
    bags = []
    for b in range(cfg.n_bags):
        neg = _draw_unlabeled(rng, unl_idx, n_pos)
        rows = np.concatenate([pos_idx, neg])
        yb = np.concatenate([np.ones(n_pos), np.zeros(len(neg))])
        clf = _lgbm(cfg, seed + b) if model == "lgbm" else _elasticnet()
        clf.fit(Xn[rows], yb)
        bags.append(clf)
        p = clf.predict_proba(Xn)[:, 1]
        in_bag_mask = np.zeros(n, dtype=bool)
        in_bag_mask[rows] = True
        oob_mask = ~in_bag_mask
        oob_sum[oob_mask] += p[oob_mask]
        oob_cnt[oob_mask] += 1
        inbag_sum[in_bag_mask] += p[in_bag_mask]
        inbag_cnt[in_bag_mask] += 1

    scores = np.where(
        oob_cnt > 0, oob_sum / np.maximum(oob_cnt, 1),
        inbag_sum / np.maximum(inbag_cnt, 1),
    )
    return scores, inbag_cnt, oob_cnt, bags
