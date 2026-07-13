# Third-party licenses and data provenance

This project redistributes a small set of metadata/input files (under
`inputs/`) that originate from third-party sources. Their licenses and
attribution are reproduced below, as required.

---

## 1. GWT_perturbseq_analysis_2025 (code + assembled inputs)

Source: https://github.com/emdann/GWT_perturbseq_analysis_2025
Upstream commit vendored here: `848d62fc2b7027f7218d6fc5f5b0c37255dc94af`

The following files in `inputs/` were taken from this repository, and the
`regulator_burden_correlation` scoring in `iei_regulator_burden_pipeline.py`
is a Python reimplementation of that repository's reference R script
(`src/8_lymphocyte_counts_LoF/Regulator_burden_correlation_GWT.R`):

  - IUIS-IEI-list-July-2024V2.csv
  - sgRNA_library_curated.csv
  - DE_stats.suppl_table.csv
  - Backman_LymphocyteCount_fullFeatures.per_gene_estimates.tsv
  - Backman_2021_86_fullFeatures.per_gene_estimates.tsv
  - shet_10bins.txt
  - gencode_v41_gname_gid_ALL_sorted
  - core_genes_reference/core_genes_*.txt

This repository is distributed under the MIT License:

```
MIT License

Copyright (c) 2025 Emma Dann

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## 2. Upstream data sources for the vendored inputs

The MIT license above covers the GWT repository. Several of the vendored
files are *data* whose ultimate provenance and terms flow from upstream
sources — cite these as well:

- **IUIS-IEI-list-July-2024V2.csv** — derived from the International Union of
  Immunological Societies (IUIS) Expert Committee classification of inborn
  errors of immunity (2024 update). Cite the IUIS IEI committee classification.

- **Backman_*_fullFeatures.per_gene_estimates.tsv** — per-gene loss-of-function
  burden effect sizes on UK Biobank quantitative traits, GeneBayes-denoised.
  Underlying burden statistics: Backman et al., "Exome sequencing and analysis
  of 454,787 UK Biobank participants," *Nature* 599, 628-634 (2021), deposited
  in the GWAS Catalog. Denoising via GeneBayes (Zeng et al.). Assembled by the
  GWT authors.

- **shet_10bins.txt** — gene-level selective-constraint (s_het) estimates.
  Cite the source constraint estimates used by the GWT study.

- **gencode_v41_gname_gid_ALL_sorted** — Ensembl gene id <-> symbol map from
  GENCODE release 41. Cite GENCODE / Ensembl.

---

## 3. Differential-expression matrix (not redistributed)

The DE effect matrix (`log_fc`) is streamed at runtime from the CZI Virtual
Cell Models bucket and is **not** redistributed by this project:

  s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad

Cite the genome-scale CD4+ T-cell Perturb-seq study and the CZI Virtual Cell
Models dataset when using it.
