# Data provenance

This repo commits only **project-derived intermediates** and **small subsets** of
large public datasets. No large raw file is redistributed. Every external source
is listed below with how to obtain it and how it feeds the pipeline.

Two reproduction tiers:

- **Tier 1 (default, offline):** reproduce `outputs/scored_full_gene_list.tsv` and
  every figure in `plots/` from the committed data — no downloads, no network.
- **Tier 2 (full recompute):** regenerate the derived intermediates from their
  original sources. Needs the Open Targets API and the S3 perturb-seq `.h5ad`.

---

## Committed — project-derived (exist nowhere else)

Small intermediates this project produces; needed for Tier-1 reproduction.

| File | Produced by |
|---|---|
| `data/derived/full_gene_list.tsv` | `src/data_build/creating_full_gene_list.py` |
| `data/derived/gwas_gene_scores.csv` | `src/gwas/gwas_gene_scores.py` |
| `data/derived/regulator_burden_scores_*.csv` | `src/regulator_burden/regulator_burden_pipeline.py` |
| `data/derived/regulator_burden_scores_wide.csv` | `src/data_build/regulator_burden_wide.py` |
| `data/derived/polarization_score.csv`, `polarization_gene_scores.csv` | `src/regulator_burden/polarization_score.py` |
| `data/derived/knn_immune_target_candidates_*.csv`, `knn_signature_pcs_*.csv` | `src/knn/` |
| `data/derived/iei_enrichment_*.csv`, `gene_overlap_summary.csv` | `src/iei_enrichment/make_enrichment_table.py` |
| `data/derived/trial_validation_top25.csv`, `trial_validation_counts.json` | `src/figures/fig_trial.py --fetch` |
| `data/derived/immune_system_drugs_grouped.csv` | `src/data_build/group_by_target_action.py` |

## Committed — reference sets & small subsets

| File | Source | Notes |
|---|---|---|
| `data/reference/IEI_gene_list.csv` | IUIS Inborn Errors of Immunity classification (Tangye et al. 2022) | 504-gene curated panel |
| `data/reference/fda_approved_target_genes.txt` | Target genes of FDA-approved immune-indicated drugs (Open Targets) | 723 symbols; the PU-model supervision labels |
| `data/reference/hgnc_gene_groups.tsv` | HGNC gene groups | frozen snapshot; defines CV gene families (exact reproducibility depends on this snapshot) |
| `data/reference/druggable_genome_gene_list.xlsx` | Finan et al. 2017, *Sci. Transl. Med.* (druggable genome) | 68 KB; committed directly |
| `data/reference/gnomad_constraint_subset.tsv` | **subset** of gnomAD v4.1.1 constraint metrics (canonical) | 278 MB → 1.7 MB; see `scripts/subset_gnomad.py` |
| `data/reference/hpa_immune_max_ntpm.tsv` | **subset** of Human Protein Atlas `rna_immune_cell.tsv` | 17 MB → 0.5 MB (per-gene max nTPM); see `scripts/subset_hpa.py` |
| `data/reference/approved_target_drugs.csv` | **subset** of the raw Open Targets drug pull | 56 MB → 0.2 MB (approved target→drug names); see `scripts/subset_approved_drugs.py` |
| `src/regulator_burden/inputs/*` | `emdann/GWT_perturbseq_analysis_2025` (MIT), Backman et al. 2021, GENCODE v41 | small vendored metadata; see `src/regulator_burden/SOURCE.md` and `THIRD_PARTY_LICENSES.md` |

## Not committed — fetched for Tier-2 only

| Source | How it's obtained | Consumed by |
|---|---|---|
| Open Targets Platform (GraphQL API, v4) | pulled by `src/data_build/build_immune_drug_dataset.py`, `src/gwas/gwas_gene_scores.py`, `src/figures/fig_trial.py --fetch` | drug list, GWAS scores, trial validation |
| CD4⁺ perturb-seq DE matrix (`GWCD4i.DE_stats.h5ad`, ~16.8 GB, public S3) | streamed via HTTP range requests (no full download); cached as a ~1.3 GB `*_f32.dat` memmap (gitignored) | `regulator_burden_pipeline.py` (`log_fc`), `knn_immune_target_score.py` (`zscore`) |
| gnomAD v4.1.1 constraint metrics (canonical), full (278 MB) | download from the gnomAD downloads page into `data/raw/`; only needed to **regenerate** `gnomad_constraint_subset.tsv` | `scripts/subset_gnomad.py` |
| HPA `rna_immune_cell.tsv`, full (17 MB) | download from proteinatlas.org into `data/raw/`; only needed to **regenerate** `hpa_immune_max_ntpm.tsv` | `scripts/subset_hpa.py` |

Run `python scripts/fetch_data.py` (or `make fetch-data`) to drive the Open Targets
pulls; the S3 `.h5ad` is streamed automatically by the scoring scripts on first run.
