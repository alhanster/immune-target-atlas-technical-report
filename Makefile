# Immune Target Atlas — reproducible pipeline
#
# Tier 1 (default, offline): reproduce the ranked gene list and all figures from
#   the committed data in data/derived + data/reference. Needs Python + R only.
# Tier 2 (full recompute): regenerate the derived intermediates from original
#   sources (Open Targets API + the S3 perturb-seq .h5ad). See DATA.md.

PY := ./.venv/bin/python
PIP := ./.venv/bin/pip
export PYTHONPATH := src

.PHONY: help setup install model scored figures test clean \
        fetch-data gene-list scores tier2

help:
	@echo "Tier 1 (offline, from committed data):"
	@echo "  make setup     - create .venv and install Python deps"
	@echo "  make model     - run the PU model -> outputs/{nominations,metrics,report,shap}"
	@echo "  make scored    - join PU scores -> outputs/scored_full_gene_list.tsv"
	@echo "  make figures   - render all 7 manuscript figures -> plots/ (needs R)"
	@echo "  make test      - run the pipeline unit tests"
	@echo "  make clean     - remove generated outputs/ and plots/"
	@echo ""
	@echo "Tier 2 (full recompute from original sources; see DATA.md):"
	@echo "  make fetch-data - stream/pull the raw inputs into data/raw/"
	@echo "  make gene-list  - rebuild data/derived/full_gene_list.tsv"
	@echo "  make tier2      - regenerate all derived intermediates"

# ---- Tier 1 -----------------------------------------------------------------
setup: install
install:
	python3 -m venv .venv
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

model:
	$(PY) -m model.run

# quick smoke run (15 bags, no ablations) — for sanity, not for publication numbers
quick:
	$(PY) -m model.run --quick --skip-ablations

scored:
	$(PY) src/model/build_scored_full_gene_list.py

figures:
	Rscript src/figures/fig_iei_enrichment.R
	Rscript src/figures/fig_regulator_burden.R
	Rscript src/figures/fig_polarization_score.R
	Rscript src/figures/fig_gwas_violin.R
	Rscript src/figures/fig_knn_candidates.R
	Rscript src/figures/fig_model.R
	Rscript src/figures/fig_trial.R

test:
	$(PY) -m pytest -q tests

clean:
	rm -f outputs/nominations.tsv outputs/metrics.json outputs/report.md \
	      outputs/top_unlabeled_shap.tsv outputs/scored_full_gene_list.tsv
	rm -f plots/*.png

# ---- Tier 2 (full recompute) ------------------------------------------------
fetch-data:
	$(PY) scripts/fetch_data.py

# rebuild the master gene table from the committed gnomAD subset + score joins
gene-list:
	$(PY) src/data_build/creating_full_gene_list.py

# regenerate every derived intermediate (requires network + the S3 .h5ad stream)
tier2: fetch-data
	$(PY) src/gwas/gwas_gene_scores.py
	$(PY) src/regulator_burden/regulator_burden_pipeline.py
	$(PY) src/regulator_burden/polarization_score.py
	$(PY) src/data_build/regulator_burden_wide.py
	$(PY) src/knn/knn_immune_target_score.py
	$(PY) src/knn/embed_signatures.py
	$(PY) src/knn/add_fda_drugs_column.py
	$(PY) src/data_build/creating_full_gene_list.py
	$(PY) src/iei_enrichment/make_enrichment_table.py
