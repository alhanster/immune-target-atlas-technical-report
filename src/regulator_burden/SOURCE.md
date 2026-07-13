# Provenance of vendored input files

All files in this directory are redistributed from third-party sources.
See `../THIRD_PARTY_LICENSES.md` for full license texts and citations.

Primary source: **emdann/GWT_perturbseq_analysis_2025** (MIT License,
Copyright (c) 2025 Emma Dann)
Vendored from upstream commit: `848d62fc2b7027f7218d6fc5f5b0c37255dc94af`

The IEI gene panel is **not** from the GWT repo — it is a user-supplied list
(`IEI_gene_list.csv`, one `Gene` column of HGNC symbols) and replaces the
IUIS-list scrape used in the original design. The pipeline no longer reads it;
IEI selection now lives in `make_regulator_burden_figure.R` Panel b, which reads
the identical copy at `../Data/Gene List/IEI_gene_list.csv`.

| file (here) | original path in GWT repo | upstream origin |
|---|---|---|
| `IEI_gene_list.csv` | — (user-supplied) | curated IEI panel |
| `sgRNA_library_curated.csv` | `metadata/` | GWT study |
| `DE_stats.suppl_table.csv` | `metadata/suppl_tables/` | GWT study |
| `Backman_LymphocyteCount_fullFeatures.per_gene_estimates.tsv` | `src/8_lymphocyte_counts_LoF/input/` | Backman 2021 (UKB) + GeneBayes |
| `Backman_2021_86_fullFeatures.per_gene_estimates.tsv` | `src/8_lymphocyte_counts_LoF/input/` | Backman 2021 (UKB) + GeneBayes |
| `shet_10bins.txt` | `src/8_lymphocyte_counts_LoF/input/` | s_het constraint estimates |
| `gencode_v41_gname_gid_ALL_sorted` | `src/8_lymphocyte_counts_LoF/input/` | GENCODE v41 |
| `core_genes_reference/core_genes_*.txt` | `src/8_lymphocyte_counts_LoF/results/` | GWT study (used for validation only) |

The large differential-expression matrix (`log_fc`, the beta term) is **not**
vendored — it is streamed at runtime from the CZI Virtual Cell Models bucket.
