"""Tier-2 data fetch orchestrator.

Tier-1 reproduction (the ranked gene list + all figures) needs NOTHING from here
— it runs entirely from the committed data in data/derived and data/reference.

This script is only for a *full recompute* of the derived intermediates from
their original sources. It handles the inputs the code can genuinely pull, and
prints guidance for the ones that need a manual download or are streamed on
demand. See DATA.md for the full provenance table.

Usage:
    python scripts/fetch_data.py            # run the Open Targets pulls
    python scripts/fetch_data.py --drugs    # only the drug-list pull
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable
RAW = REPO / "data" / "raw"


def run(desc: str, argv: list[str]) -> None:
    print(f"\n=== {desc} ===")
    subprocess.run(argv, check=True, cwd=REPO, env={"PYTHONPATH": "src", **_env()})


def _env() -> dict:
    import os
    return dict(os.environ)


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    only_drugs = "--drugs" in sys.argv

    # 1. Open Targets GraphQL pulls (no manual download — the scripts hit the API).
    run("Open Targets: immune drug dataset -> data/raw/immune_system_drugs.csv",
        [PY, "src/data_build/build_immune_drug_dataset.py"])
    run("Open Targets: GWAS credible-set scores -> data/derived/gwas_gene_scores.csv",
        [PY, "src/gwas/gwas_gene_scores.py"])

    if only_drugs:
        return

    # 2. Inputs that need a manual download (committed subsets already cover the
    #    Tier-1 path; only needed to regenerate the subsets themselves).
    print("\n=== Manual downloads (only needed to regenerate committed subsets) ===")
    print("See DATA.md for exact sources. Place files under data/raw/, then run the")
    print("matching scripts/subset_*.py to refresh data/reference/:")
    print("  - gnomAD v4.1.1 constraint metrics (canonical) -> scripts/subset_gnomad.py")
    print("  - HPA rna_immune_cell.tsv                       -> scripts/subset_hpa.py")
    print("  - raw immune_system_drugs.csv (from step 1)     -> scripts/subset_approved_drugs.py")

    # 3. Streamed on demand — no explicit download.
    print("\n=== Streamed on demand (no download step) ===")
    print("The 16.8 GB CD4+ perturb-seq .h5ad is read directly from its public S3")
    print("bucket via HTTP range requests by:")
    print("  - src/regulator_burden/regulator_burden_pipeline.py  (log_fc layer)")
    print("  - src/knn/knn_immune_target_score.py                 (zscore layer)")
    print("Each caches a ~1.3 GB float32 memmap under its inputs/ or data/raw/ (gitignored).")

    print("\nDone. For the full recompute chain, run `make tier2`.")


if __name__ == "__main__":
    main()
