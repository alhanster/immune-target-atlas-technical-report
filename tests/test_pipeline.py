"""Unit tests for the PU pipeline. Run: make test  (or python -m pytest -q tests)."""
import numpy as np
import pandas as pd
import pytest

from model import config, data, labels, features, cv, evaluate
from model.config import Config


@pytest.fixture(scope="module")
def loaded():
    df = data.load_gene_table()
    fda = data.load_fda_targets()
    df = labels.build_labels(df, fda)
    family, info = data.load_family_groups(df["gene"].tolist())
    return df, fda, family, info


def test_no_rows_dropped(loaded):
    df, *_ = loaded
    assert len(df) == 18692


def test_label_sets(loaded):
    df, fda, *_ = loaded
    lab = labels.label_summary(df, fda)
    # set b (immune) must be a strict subset of set a (all)
    assert lab["positives_immune"] <= lab["positives_all"]
    assert lab["positives_all"] == int(df["label_all"].sum())
    # every immune positive is an approved target with immune evidence
    imm = df[df["label_immune"] == 1]
    assert imm["label_all"].all()
    assert ((imm["gwas_score"] > 0) | (imm["IEI"] == 1)).all()


def test_gwas_provenance_flag():
    # provenance is confirmed from Part 3 source; guardrail banner must reflect it
    assert labels.GWAS_PROVENANCE_CONFIRMED is True
    assert "gwas_credible_sets" in labels.GWAS_PROVENANCE_NOTE


def test_features_keep_all_rows_and_indicators(loaded):
    df, _, family, _ = loaded
    cfg = Config()
    X, groups, presence = features.build_features(df, cfg, family=family)
    assert len(X) == len(df)
    # presence indicators exist and are 0/1
    for c in ["gwas_present", "perturbseq_measured", "polarization_measured"]:
        assert set(presence[c].unique()) <= {0, 1}
    # impute0 mode -> no NaN in assay columns
    assert not X["rest_signed_log10p"].isna().any()
    # constraint columns keep native NaN (MAR)
    assert X["lof.oe_ci.upper"].isna().any()


def test_native_mode_leaves_assay_nan(loaded):
    df, _, family, _ = loaded
    cfg = Config(assay_na_mode="native")
    X, groups, _ = features.build_features(df, cfg, family=family)
    assert X["rest_signed_log10p"].isna().any()
    assert "perturbseq_measured" not in groups["perturbseq"]


def test_grouped_folds_no_family_leak(loaded):
    df, _, family, _ = loaded
    y = df["label_immune"].to_numpy()
    folds = cv.make_folds(family, y, n_folds=5, seed=0)
    groups = family.to_numpy()
    for tr, ev in folds:
        assert not (set(groups[tr]) & set(groups[ev]))
        assert y[tr].sum() > 0 and y[ev].sum() > 0  # positives in every split


def test_pu_metrics_perfect_ranking():
    y = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    scores = np.array([9, 8, 1, 1, 1, 1, 1, 1, 1, 1], dtype=float)
    m = evaluate.pu_metrics(y, scores, top_k_list=(2,))
    assert m["recall@2"] == 1.0
    assert m["enrichment@2"] == pytest.approx(1.0 / m["base_rate"])
