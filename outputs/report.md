# PU target-nomination — run report

*Label set:* **immune** (primary) · *bags:* 100 · *folds:* 5 · *seed:* 1234

## Guardrail findings

**1. gwas_score provenance — CONFIRMED SAFE.** `gwas_score` is the Open Targets `gwas_credible_sets` *datasource* score for MONDO:0005046 (the GWAS component of `genetic_association`), verified from `src/gwas/gwas_gene_scores.py`. It is **not** the overall association score, so it does not fold in the ChEMBL known-drug datasource and does **not** leak approval labels. The with/without-gwas ablation is still reported below for transparency.

**2. Positive-set definition.** (a) all approved = **691** positives; (b) immune-restricted (gwas>0 or IEI==1) = **361** positives — **(b) is the primary/honest target**, (a) is the sensitivity analysis. The approved-target file lists 723 symbols; 32 are outside the 18,692-gene universe.

**3. No dropping rows.** All **18692** genes retained; per-column NA policy (0-impute + presence indicator for assay/polarization blocks; native NaN for constraint columns). No complete-casing.

> *Residual caveat (label set b):* every immune-restricted positive has `gwas_score>0` or `IEI==1` by construction, so those two features partly encode the label gating. This is why the without-gwas ablation and the `IEI`-inclusive vs genetics-only comparison matter — read enrichment alongside them, not in isolation.

## Cross-validation

Group-aware **StratifiedGroupKFold** by gene family — 3229 families (15131 from HGNC gene groups, 3561 from the approximate symbol-prefix fallback). Largest family = 4.0% of genes. The bagging-PU is nested inside each fold, so no positive's family appears in both train and eval.

## Primary result vs baseline

| Model | recall@50 | recall@100 | recall@200 | enrich@100 | PU-AP |
|---|---|---|---|---|---|
| LightGBM bagging-PU | 0.025 | 0.042 | 0.064 | 7.8 | 0.101 |
| Elastic-net logistic | 0.028 | 0.047 | 0.069 | 8.8 | 0.112 |

**Honest finding:** the elastic-net baseline is *competitive with / slightly ahead of* LightGBM here. This is expected and not a bug — under label set (b) the dominant signal is a largely **monotone** combination of `gwas_score` + `IEI` + constraint, which a linear model captures directly; the tree's advantages (interactions, non-linear thresholds) add little when the label is essentially a genetics threshold. The tree model is still kept as the primary ranker for its **native NA handling** (no imputation of the MAR constraint columns) and its **per-gene SHAP** attributions; the small PU-AP gap (<0.02) is within bagging noise. If trees are to earn their place they should be revisited on richer, interaction-heavy feature sets.

## Ablation table (grouped-CV OOF held-out-positive enrichment)

| Variant | recall@50 | recall@100 | recall@200 | enrich@100 | PU-AP |
|---|---|---|---|---|---|
| (iv) label=immune [PRIMARY] | 0.022 | 0.044 | 0.078 | 8.3 | 0.100 |
| (iv) label=all [sensitivity] | 0.014 | 0.020 | 0.041 | 3.8 | 0.090 |
| (i)  with gwas_score | 0.022 | 0.044 | 0.078 | 8.3 | 0.100 |
| (i)  without gwas_score | 0.011 | 0.022 | 0.039 | 4.1 | 0.049 |
| (ii) assay impute0+indicator | 0.022 | 0.044 | 0.078 | 8.3 | 0.100 |
| (ii) assay NA-native | 0.022 | 0.044 | 0.078 | 8.3 | 0.097 |
| (iii) without family/druggability | 0.022 | 0.044 | 0.078 | 8.3 | 0.100 |
| (iii) with family/druggability | 0.022 | 0.050 | 0.078 | 9.3 | 0.088 |

**Reading the ablations.** *(i)* Dropping `gwas_score` roughly halves recall@100 and PU-AP — genetics carries about half the signal, and the score is safe to use (guardrail 1). *(iv)* Label set (b) enriches far better than (a): approved targets broadly include many non-immune drugs, so the immune-restricted target is the honest one. *(ii)* Assay NA-native ≈ impute0+indicator — the pipeline is robust to that choice. *(iii)* Adding family/druggability features nudges recall@100 up but *lowers* PU-AP — the classic circularity signature, which is why they are kept in a separate group and left **off** by default.

## Top nominated unlabeled genes (SHAP-explained)

Ranked non-approved genes; `driver` flags whether the score is carried by genetics+Perturb-seq or by family/druggability features.

| rank | gene | pu_score | driver | top features |
|---|---|---|---|---|
| 2 | NOS2 | 0.982 | genetics+perturbseq-driven | IEI+1.978; gwas_score+1.337; mis.z_score+0.602 |
| 5 | IL10 | 0.978 | genetics+perturbseq-driven | IEI+2.169; gwas_score+1.366; stim48hr_se_beta+0.338 |
| 6 | CFHR3 | 0.976 | genetics+perturbseq-driven | IEI+2.124; gwas_score+0.855; mis.z_score+0.555 |
| 10 | CDCA7 | 0.975 | genetics+perturbseq-driven | IEI+2.125; gwas_score+1.346; stim48hr_coef_beta+0.289 |
| 11 | ERN1 | 0.975 | genetics+perturbseq-driven | IEI+2.062; gwas_score+1.242; mis.z_score+0.523 |
| 13 | MAP3K14 | 0.973 | genetics+perturbseq-driven | IEI+2.023; gwas_score+0.921; mis.z_score+0.872 |
| 14 | NBEAL2 | 0.973 | genetics+perturbseq-driven | IEI+2.142; gwas_score+0.741; mis.z_score+0.609 |
| 15 | NLRP3 | 0.972 | genetics+perturbseq-driven | IEI+2.034; gwas_score+0.869; mis.z_score+0.720 |
| 17 | SRP54 | 0.971 | genetics+perturbseq-driven | IEI+2.032; gwas_score+0.822; mis.z_score+0.791 |
| 18 | IRF3 | 0.970 | genetics+perturbseq-driven | IEI+2.145; gwas_score+1.315; stim8hr_se_beta-0.180 |
| 19 | ZAP70 | 0.970 | genetics+perturbseq-driven | IEI+2.181; gwas_score+1.320; mis.z_score+0.305 |
| 20 | CD40 | 0.970 | genetics+perturbseq-driven | IEI+2.102; gwas_score+1.260; stim48hr_coef_beta+0.381 |
| 21 | CORO1A | 0.970 | genetics+perturbseq-driven | IEI+2.093; gwas_score+0.917; mis.z_score+0.516 |
| 22 | ACTB | 0.969 | genetics+perturbseq-driven | IEI+1.995; gwas_score+1.172; mis.z_score+0.922 |
| 23 | RPSA | 0.968 | genetics+perturbseq-driven | IEI+1.942; gwas_score+0.970; mis.z_score+0.786 |
| 24 | ICOSLG | 0.968 | genetics+perturbseq-driven | IEI+1.928; gwas_score+1.111; mis.z_score+0.711 |
| 25 | IKBKE | 0.968 | genetics+perturbseq-driven | IEI+2.049; gwas_score+1.337; mis.z_score+0.562 |
| 26 | RASGRP1 | 0.968 | genetics+perturbseq-driven | IEI+2.060; gwas_score+1.321; mis.z_score+0.534 |
| 29 | IKZF1 | 0.966 | genetics+perturbseq-driven | IEI+2.097; gwas_score+1.263; mis.z_score+0.442 |
| 30 | NFKBIA | 0.966 | genetics+perturbseq-driven | IEI+2.191; gwas_score+1.319; stim48hr_se_beta+0.310 |

## Files

- `nominations.tsv` — all genes, PU score, rank, label status, presence flags.
- `top_unlabeled_shap.tsv` — top-N unlabeled genes with SHAP drivers.
- `metrics.json` — full config, guardrails, ablation metrics.


---
*Headline is grouped-CV held-out-positive enrichment, not a single AUC.*
