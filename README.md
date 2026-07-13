# Immune Target Atlas

A four-layer human-genetics and functional-genomics framework that prioritizes
novel drug targets for immune and autoimmune disease. It integrates four
orthogonal evidence layers per gene and fuses them into a single genome-wide
ranking of all 18,692 protein-coding genes:

1. **Monogenic ground truth** — inborn errors of immunity (IUIS IEI panel).
2. **Constraint** — gnomAD evolutionary/population constraint (missense-z, LOEUF).
3. **Complex-trait genetics** — Open Targets GWAS credible-set score for immune
   system disorder (MONDO:0005046).
4. **Causal functional genomics** — genome-scale CRISPRi perturb-seq in primary
   human CD4⁺ T cells (regulator-burden + Th1/Th2 polarization + nearest-target kNN).

A **positive–unlabeled (PU) learning** model treats FDA-approved immune-drug
targets as positives and ranks every gene by immune-target potential, with
group-aware cross-validation (gene-family holdout) and per-gene SHAP attribution.

## Repository layout

```
src/
  shared/            gene_name_utils.py — HGNC symbol harmonization (shared)
  data_build/        build the master gene table + drug-list ingestion
  iei_enrichment/    Part 1 — IEI vs FDA-target Fisher enrichment
  gwas/              Part 3 — Open Targets GWAS credible-set scores
  regulator_burden/  Part 2 — regulator-burden + Th1/Th2 polarization (perturb-seq)
  knn/               Part 4 — nearest-FDA-target cosine similarity in signature space
  model/             the PU-learning ranker (run: python -m model.run)
  figures/           R (+ Python) renderers for the 7 manuscript figures
data/
  reference/         curated sets + small committed subsets of large sources
  derived/           project-produced intermediates (committed)
  raw/               gitignored; Tier-2 downloads land here
outputs/             nominations.tsv, scored_full_gene_list.tsv, metrics.json, report.md, SHAP
plots/               regenerated manuscript figures (PNG)
scripts/             fetch_data.py + subset_*.py helpers
docs/                methods notes; docs/paper/ reserved for the manuscript
tests/               pipeline unit tests
```

## Quickstart (Tier 1 — offline, from committed data)

Reproduces the ranked gene list and all figures with no downloads.

```bash
make setup      # create .venv, install Python deps (see requirements.txt)
make model      # PU model -> outputs/{nominations,metrics,report,top_unlabeled_shap}
make scored     # -> outputs/scored_full_gene_list.tsv
make test       # pipeline unit tests

# figures need R (see R-requirements.txt):
#   install.packages(c("ggplot2","patchwork","jsonlite","ragg","uwot","scales"))
make figures    # -> plots/*.png
```

The model is deterministic (seed 1234): `outputs/scored_full_gene_list.tsv` and the
figures reproduce exactly, given the same LightGBM version.

### Key results (from `outputs/`)

- PU model recovers held-out known targets at **7.8× enrichment** in its top 100.
- Removing the GWAS layer roughly halves enrichment (**8.3× → 4.1×**) — human
  genetics carries about half the signal.
- Highlighted novel candidates: **MAP3K14 (NIK), RORC (RORγt), MALT1** — all on the
  antigen-receptor→NF-κB activation axis, none yet an approved immune-drug target.

## Full recompute (Tier 2)

Regenerates the derived intermediates from their original sources (Open Targets
API + the 16.8 GB S3 perturb-seq `.h5ad`, streamed on demand). See **[DATA.md](DATA.md)**
for the provenance table.

```bash
make fetch-data   # Open Targets pulls into data/raw/ (S3 .h5ad streams on first use)
make tier2        # regenerate every data/derived/ intermediate
```

## Data & licensing

- Data provenance and the commit-vs-fetch policy: **[DATA.md](DATA.md)**.
- Code is MIT licensed (**[LICENSE](LICENSE)**). Part-2 regulator-burden code is
  derived from `emdann/GWT_perturbseq_analysis_2025` (MIT); see
  `src/regulator_burden/SOURCE.md` and `THIRD_PARTY_LICENSES.md`.
- All primary data derive from public resources: IUIS IEI, gnomAD, Open Targets,
  Human Protein Atlas, Finan et al. 2017, and the CD4⁺ perturb-seq dataset of
  Zhu et al. 2025.
