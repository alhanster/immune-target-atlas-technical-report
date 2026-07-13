"""End-to-end CLI: build features, run grouped-CV enrichment, ablations, and
produce the final PU nominations with SHAP explanations.

    python -m model.run                 # full run, default (honest) config
    python -m model.run --quick         # fewer bags, for a fast smoke run
    python -m model.run --help
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import platform
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

from . import config, data, labels, features, evaluate, model as model_mod, shap_explain
from .config import Config


def _y_for(df: pd.DataFrame, label_set: str) -> np.ndarray:
    col = "label_immune" if label_set == "immune" else "label_all"
    return df[col].to_numpy()


def run_cv_variant(df, family, cfg: Config, model: str = "lgbm") -> dict:
    """Grouped-CV OOF PU metrics for one config variant."""
    X, groups, _ = features.build_features(df, cfg, family=family)
    cols = features.active_feature_columns(groups, cfg.use_druggability)
    Xc = X[cols]
    y = _y_for(df, cfg.label_set)
    oof = evaluate.cross_val_oof_scores(Xc, y, family, cfg, model=model)
    m = evaluate.pu_metrics(y, oof, cfg.top_k_list)
    m["model"] = model
    m["n_features"] = len(cols)
    return m


def build_ablation_table(df, family, base: Config, ablation_bags: int) -> list[dict]:
    """Run the four required ablations, each as grouped-CV OOF enrichment."""
    def cfg(**kw):
        return dataclasses.replace(base, n_bags=ablation_bags, **kw)

    rows = []
    variants = [
        ("(iv) label=immune [PRIMARY]", cfg(label_set="immune")),
        ("(iv) label=all [sensitivity]", cfg(label_set="all")),
        ("(i)  with gwas_score", cfg(use_gwas=True)),
        ("(i)  without gwas_score", cfg(use_gwas=False)),
        ("(ii) assay impute0+indicator", cfg(assay_na_mode="impute0_indicator")),
        ("(ii) assay NA-native", cfg(assay_na_mode="native")),
        ("(iii) without family/druggability", cfg(use_druggability=False)),
        ("(iii) with family/druggability", cfg(use_druggability=True)),
    ]
    for name, c in variants:
        t0 = time.time()
        m = run_cv_variant(df, family, c, model="lgbm")
        m["variant"] = name
        m["seconds"] = round(time.time() - t0, 1)
        rows.append(m)
        print(f"  [{name:38s}] recall@100={m.get('recall@100', float('nan')):.3f} "
              f"enrich@100={m.get('enrichment@100', float('nan')):.1f} "
              f"PU-AP={m['pu_average_precision']:.3f}  ({m['seconds']}s)")
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description="PU target-nomination pipeline")
    ap.add_argument("--label-set", choices=["immune", "all"], default="immune")
    ap.add_argument("--no-gwas", action="store_true", help="drop gwas_score")
    ap.add_argument("--assay-na", choices=["impute0_indicator", "native"],
                    default="impute0_indicator")
    ap.add_argument("--druggability", action="store_true",
                    help="add ablatable family/druggability features to core model")
    ap.add_argument("--drop-coef-se", action="store_true")
    ap.add_argument("--n-bags", type=int, default=100)
    ap.add_argument("--ablation-bags", type=int, default=40,
                    help="bags per ablation variant (speed vs stability)")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--quick", action="store_true",
                    help="fast smoke run: 15 bags, 8 ablation bags")
    ap.add_argument("--skip-ablations", action="store_true")
    args = ap.parse_args(argv)

    n_bags = 15 if args.quick else args.n_bags
    ablation_bags = 8 if args.quick else args.ablation_bags

    cfg = Config(
        use_gwas=not args.no_gwas,
        label_set=args.label_set,
        assay_na_mode=args.assay_na,
        drop_coef_se=args.drop_coef_se,
        use_druggability=args.druggability,
        n_bags=n_bags,
        n_folds=args.folds,
        seed=args.seed,
        top_n_report=args.top_n,
    )
    cfg.out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load -----------------------------------------------------------
    print("Loading data ...")
    df = data.load_gene_table()
    fda = data.load_fda_targets()
    df = labels.build_labels(df, fda)
    family, fam_info = data.load_family_groups(df["gene"].tolist())
    assert len(df) == 18692, f"expected 18,692 genes, got {len(df)} (rows must not drop)"

    lab = labels.label_summary(df, fda)

    # ---- Guardrail banners ---------------------------------------------
    print("\n" + "=" * 78)
    print(labels.gwas_provenance_warning())
    print("-" * 78)
    print("[GUARDRAIL 2 — positive-set definition]")
    print(f"  (a) all approved targets in gene list : {lab['positives_all']} positives")
    print(f"  (b) immune-restricted (gwas>0 or IEI) : {lab['positives_immune']} positives  <- PRIMARY")
    print(f"  FDA file lists {lab['n_fda_file']}; {lab['n_fda_missing_from_list']} not in the 18,692-gene universe.")
    print("-" * 78)
    print("[GUARDRAIL 3 — no dropping rows]")
    print(f"  Kept all {lab['n_genes']} genes. NA handled per-column (see README); "
          "no complete-casing.")
    print("-" * 78)
    print(f"[CV grouping] {fam_info['n_families']} families | "
          f"{fam_info['n_from_hgnc']} HGNC, {fam_info['n_from_prefix_fallback']} prefix-fallback "
          f"(approximate) | largest family = {fam_info['largest_family_frac']*100:.1f}% of genes")
    print("=" * 78 + "\n")

    # ---- Primary grouped-CV result + baseline ---------------------------
    print(f"Primary grouped-CV enrichment  (label={cfg.label_set}, bags={cfg.n_bags}) ...")
    primary = run_cv_variant(df, family, cfg, model="lgbm")
    print(f"  LightGBM   : recall@100={primary['recall@100']:.3f} "
          f"enrich@100={primary['enrichment@100']:.1f} PU-AP={primary['pu_average_precision']:.3f}")
    print("Elastic-net logistic baseline (must underperform to justify trees) ...")
    base_cfg = dataclasses.replace(cfg, n_bags=min(cfg.n_bags, 40))
    baseline = run_cv_variant(df, family, base_cfg, model="elasticnet")
    print(f"  ElasticNet : recall@100={baseline['recall@100']:.3f} "
          f"enrich@100={baseline['enrichment@100']:.1f} PU-AP={baseline['pu_average_precision']:.3f}")

    # ---- Ablation table -------------------------------------------------
    ablation = []
    if not args.skip_ablations:
        print(f"\nAblation table (grouped-CV OOF, {ablation_bags} bags each):")
        ablation = build_ablation_table(df, family, cfg, ablation_bags)

    # ---- Final full-data nominations + SHAP -----------------------------
    print("\nFinal bagging-PU on all genes (OOB-averaged scores) ...")
    X, groups, presence = features.build_features(df, cfg, family=family)
    cols = features.active_feature_columns(groups, cfg.use_druggability)
    Xc = X[cols]
    y = _y_for(df, cfg.label_set)
    scores, inbag, oob, bags = model_mod.bagging_pu_full(Xc, y, cfg, cfg.seed, "lgbm")

    nom = pd.DataFrame({
        "gene": df["gene"].values,
        "gene_id": df["gene_id"].values,
        "pu_score": scores,
        "label_all": df["label_all"].values,
        "label_immune": df["label_immune"].values,
        "label_status": np.where(y == 1, "positive", "unlabeled"),
        "family": family.values,
        "gwas_present": presence["gwas_present"].values,
        "perturbseq_measured": presence["perturbseq_measured"].values,
        "polarization_measured": presence["polarization_measured"].values,
        "oob_models": oob.astype(int),
    })
    nom = nom.sort_values("pu_score", ascending=False).reset_index(drop=True)
    nom.insert(0, "rank", np.arange(1, len(nom) + 1))
    nom_path = cfg.out_dir / "nominations.tsv"
    nom.to_csv(nom_path, sep="\t", index=False, float_format="%.6f")
    print(f"  wrote {nom_path} ({len(nom)} genes)")

    # SHAP for top-N UNLABELED genes
    unlabeled_rank = nom[nom["label_status"] == "unlabeled"].head(cfg.top_n_report)
    # map back to positional index in the unsorted X
    gene_to_pos = {g: i for i, g in enumerate(df["gene"].values)}
    row_pos = np.array([gene_to_pos[g] for g in unlabeled_rank["gene"]])
    print(f"Computing SHAP for top-{len(row_pos)} unlabeled genes ...")
    shap_df = shap_explain.explain_top_genes(bags, X, cols, row_pos)

    recs = []
    for (gene, rnk, sc), (_, srow) in zip(
        unlabeled_rank[["gene", "rank", "pu_score"]].itertuples(index=False, name=None),
        shap_df.iterrows(),
    ):
        recs.append({
            "rank": int(rnk),
            "gene": gene,
            "pu_score": round(float(sc), 6),
            "driver": shap_explain.driver_tag(srow),
            "top_features": shap_explain.top_feature_string(srow),
            "gwas_present": int(presence["gwas_present"].iloc[gene_to_pos[gene]]),
            "perturbseq_measured": int(presence["perturbseq_measured"].iloc[gene_to_pos[gene]]),
        })
    shap_out = pd.DataFrame(recs)
    shap_path = cfg.out_dir / "top_unlabeled_shap.tsv"
    shap_out.to_csv(shap_path, sep="\t", index=False)
    print(f"  wrote {shap_path}")

    # ---- metrics.json ---------------------------------------------------
    metrics = {
        "config": cfg.to_dict(),
        "environment": {
            "python": platform.python_version(),
            "seed": cfg.seed,
        },
        "guardrails": {
            "gwas_provenance_confirmed": labels.GWAS_PROVENANCE_CONFIRMED,
            "gwas_provenance_note": labels.GWAS_PROVENANCE_NOTE,
            "label_summary": lab,
            "rows_kept": len(df),
            "cv_grouping": fam_info,
        },
        "primary": primary,
        "baseline_elasticnet": baseline,
        "ablations": ablation,
    }
    metrics_path = cfg.out_dir / "metrics.json"
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    print(f"  wrote {metrics_path}")

    _write_report(cfg, lab, fam_info, primary, baseline, ablation, shap_out, nom)
    print("\nDone.")


def _write_report(cfg, lab, fam_info, primary, baseline, ablation, shap_out, nom):
    from .report import render_report
    path = cfg.out_dir / "report.md"
    path.write_text(render_report(cfg, lab, fam_info, primary, baseline, ablation, shap_out, nom))
    print(f"  wrote {path}")


if __name__ == "__main__":
    main()
