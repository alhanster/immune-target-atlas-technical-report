"""Central configuration for the PU target-nomination pipeline.

Every knob that changes a modelling decision lives here so a run is fully
described by one Config object (which is serialised into metrics.json).
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Literal

# Repo layout ---------------------------------------------------------------
# config.py lives at <repo>/src/model/config.py -> parents[2] is the repo root.
ROOT = Path(__file__).resolve().parents[2]
GENE_LIST_TSV = ROOT / "data" / "derived" / "full_gene_list.tsv"
APPROVED_TXT = ROOT / "data" / "reference" / "approved_target_genes.txt"
HGNC_GROUPS_TSV = ROOT / "data" / "reference" / "hgnc_gene_groups.tsv"
OUT_DIR = ROOT / "outputs"

# The disease used to build gwas_score (Open Targets MONDO:0005046).
DISEASE_EFO = "MONDO_0005046"

# Column groups -------------------------------------------------------------
PERTURBSEQ_COLS = [
    "rest_coef_beta", "rest_se_beta", "rest_signed_log10p",
    "stim8hr_coef_beta", "stim8hr_se_beta", "stim8hr_signed_log10p",
    "stim48hr_coef_beta", "stim48hr_se_beta", "stim48hr_signed_log10p",
]
SIGNED_LOG10P_COLS = ["rest_signed_log10p", "stim8hr_signed_log10p", "stim48hr_signed_log10p"]
COEF_SE_COLS = [
    "rest_coef_beta", "rest_se_beta",
    "stim8hr_coef_beta", "stim8hr_se_beta",
    "stim48hr_coef_beta", "stim48hr_se_beta",
]
CONSTRAINT_COLS = ["lof.oe_ci.upper", "mis.z_score"]


@dataclass
class Config:
    # --- Guardrail-driven switches ---------------------------------------
    use_gwas: bool = True
    # 'immune' (set b: approved AND (gwas>0 OR IEI==1)) is the primary/honest
    # label; 'all' (set a: every approved target) is the sensitivity analysis.
    label_set: Literal["immune", "all"] = "immune"

    # --- NA / imputation policy ------------------------------------------
    # 'impute0_indicator': assay blocks -> fill 0 + presence indicator.
    # 'native':            leave assay NAs for LightGBM's native handling.
    assay_na_mode: Literal["impute0_indicator", "native"] = "impute0_indicator"
    polarization_impute: Literal["zero", "median"] = "zero"

    # --- Feature engineering ---------------------------------------------
    drop_coef_se: bool = False           # drop redundant coef/se, keep signed_log10p
    use_druggability: bool = False       # add ablatable druggability/family group

    # --- Model -----------------------------------------------------------
    n_bags: int = 100                    # bagging-PU base models
    num_leaves: int = 31
    min_child_samples: int = 50
    learning_rate: float = 0.05
    n_estimators: int = 200
    subsample: float = 0.8
    colsample_bytree: float = 0.8

    # --- Cross-validation -------------------------------------------------
    n_folds: int = 5
    top_k_list: tuple[int, ...] = (50, 100, 200, 500)

    # --- Reproducibility / output ----------------------------------------
    seed: int = 1234
    top_n_report: int = 50               # #unlabeled genes to SHAP-explain
    out_dir: Path = field(default=OUT_DIR)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["out_dir"] = str(self.out_dir)
        d["top_k_list"] = list(self.top_k_list)
        return d
