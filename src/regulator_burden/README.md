# IEI regulator–burden core-gene scoring

Nominates **inborn-errors-of-immunity (IEI) genes** as candidate *core genes*
for lymphocyte count, by integrating a CD4⁺ T-cell genome-scale Perturb-seq
differential-expression matrix with UK Biobank loss-of-function (LoF) burden
estimates.

This is a lightweight, self-contained reimplementation of the regulator–burden
correlation from the genome-wide T-cell Perturb-seq study. The scoring is a
Python port of the authors' reference R script
(`src/8_lymphocyte_counts_LoF/Regulator_burden_correlation_GWT.R` in the
[analysis repo](https://github.com/emdann/GWT_perturbseq_analysis_2025)).

## Method

For each candidate core gene *j*, in each culture condition, fit across all
regulators *x* (with the self-term *x = j* masked):

```
scale(γ_x)  ~  a · scale(β_{x→j})  +  b · S_het,x  +  intercept
```

| term | meaning | source |
|------|---------|--------|
| `β_{x→j}` | log₂FC of gene *j* after CRISPRi knockdown of regulator *x* | Perturb-seq DE matrix |
| `γ_x` | GeneBayes-denoised LoF **burden effect on lymphocyte count** (not the raw count, not perturb-seq) | UKB / Backman |
| `S_het,x` | selective-constraint nuisance covariate | published estimates |

The per-gene score for *j* is the coefficient **`a` (`coef_beta`)**: its **sign**
gives direction (positive → up-regulation of *j* promotes the trait) and its
**p-value** its significance.

## Quick start

```bash
python -m pip install -r requirements.txt

python regulator_burden_pipeline.py --self-test   # seconds, no data
python regulator_burden_pipeline.py               # full genome-wide scoring
```

Scoring is genome-wide — no gene panel is involved. Selecting IEI genes (or any
other set) is done by the consumer: the figure script's Panel b reads
`../Data/Gene List/IEI_gene_list.csv` and keeps those genes. See
**Selecting a gene set** below.

The small metadata inputs (γ, S_het, gene map, validation lists) are
**vendored** in `inputs/` — no clone needed. Their provenance and licenses are
in `inputs/SOURCE.md` and `THIRD_PARTY_LICENSES.md`. Only the large DE matrix
(β) is fetched at runtime.

## Runtime

The only slow step is a **one-time** stream of the DE matrix (β), which exists
only in a single 16.8 GB `.h5ad`. The script reads that file's HDF5 header over
HTTP range requests, streams **only** the contiguous `log_fc` byte range
(~2.8 GB) and caches it to disk as a float32 memmap (`logfc_f32.dat`, ~1.4 GB).
The full file is never downloaded or held in RAM.

| command | first run | cached re-run |
|---------|-----------|---------------|
| `--self-test` | seconds (no network) | — |
| default (genome-wide) | ~3 min | ~40 s |

Delete `logfc_f32.dat` or pass `--refresh` to re-stream. Needs ~2 GB free disk
and ~2 GB RAM.

## Data sources

1. **Metadata inputs** (vendored in `inputs/`) — γ (LoF burden), S_het, the
   Ensembl↔symbol map, and the authors' own core-gene lists (validation) are
   used by the pipeline. The sgRNA library, `DE_stats.suppl_table.csv`, and a
   copy of `IEI_gene_list.csv` remain vendored but are no longer read by the
   pipeline. Everything here is redistributed from the MIT-licensed
   `emdann/GWT_perturbseq_analysis_2025` and its upstream data sources; see
   `inputs/SOURCE.md` and `THIRD_PARTY_LICENSES.md`.
2. **DE matrix (β)** — streamed at runtime, CZI Virtual Cell Models bucket,
   public anonymous read (not redistributed here):
   `s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad`

## Outputs

| file | contents |
|------|----------|
| `regulator_burden_scores_{Rest,Stim8hr,Stim48hr}.csv` | genome-wide per-gene scores |
| `regulator_burden_scores_all.csv` | **main deliverable** — genome-wide scores, all conditions stacked (coef, p, rank, direction) |

The Python pipeline writes CSVs only. The manuscript figure
(`regulator_burden_manuscript.png`, signed-QQ + IEI core-gene bar chart) is
drawn by `make_regulator_burden_figure.R`, which reads the CSVs above and
selects the IEI genes for Panel b.

## Validation

In the default (genome-wide) run, the script compares its top-50 gene lists per
condition/direction against the authors' published core-gene files under
`src/8_lymphocyte_counts_LoF/results/`. Small differences (a few genes)
arise from within-subset standardization and gene-ID mapping; the top hits and
the condition-specificity of signal reproduce.

## Selecting a gene set

The pipeline scores every gene genome-wide and has no notion of a gene panel.
Any gene-set selection happens downstream, in `make_regulator_burden_figure.R`
Panel b, which reads `../Data/Gene List/IEI_gene_list.csv` (a single `Gene`
column of HGNC symbols) and keeps only those genes. To use a different set,
point that read at another CSV of the same format, or filter
`regulator_burden_scores_all.csv` by `gene` yourself.

Legacy aliases are mapped onto the DE-matrix (GENCODE v41) naming via the
`ALIAS` dict in the script; add entries there if your symbols don't match.

## Swapping the trait

γ is the only trait-specific input. To score a different trait, supply that
trait's GeneBayes-denoised per-gene burden file (same format:
`ensg`, `post_mean`, …) in place of the lymphocyte-count file. The DE matrix,
S_het, and the entire model stay unchanged. Note that the Perturb-seq cell type
(CD4⁺ T cell) must plausibly be relevant to the trait, or the correlation is
uninterpretable.

## Attribution & citation

This project reuses code and data from third parties. Full license texts and
per-file provenance: **`THIRD_PARTY_LICENSES.md`** and **`inputs/SOURCE.md`**.

- **Method & reference implementation** — `emdann/GWT_perturbseq_analysis_2025`
  (MIT License, © 2025 Emma Dann). The regulator–burden scoring here is a
  Python port of that repo's `Regulator_burden_correlation_GWT.R`. Vendored
  inputs and validation lists come from commit `848d62fc`.
- **LoF burden data** — Backman et al., *Exome sequencing and analysis of
  454,787 UK Biobank participants*, Nature 599, 628–634 (2021); GWAS Catalog.
  GeneBayes-denoised (Zeng et al.).
- **IEI gene panel** — IUIS Expert Committee classification of inborn errors of
  immunity (2024 update).
- **Gene annotation** — GENCODE release 41 / Ensembl.
- **DE matrix** — genome-scale CD4⁺ T-cell Perturb-seq study; CZI Virtual Cell
  Models dataset.

Please cite all of the above when using these scores. This repository itself is
an independent reproduction for IEI therapeutic-target triage.

*Note: this documents the applicable licenses and standard attribution
practice; it is not legal advice. If in doubt, consult your institution's
licensing office.*
