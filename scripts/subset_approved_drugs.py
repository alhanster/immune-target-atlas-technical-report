"""Subset the raw Open Targets immune-drug pull to a tiny target->approved-drugs
table, so the repo ships a small committed file instead of the 56 MB raw pull.

The full ``immune_system_drugs.csv`` (~324k rows) is only consumed by
``src/knn/add_approved_drugs_column.py``, which needs, per target gene, the set of
approved drug names (``maxClinicalStage == "APPROVAL"``). This keeps exactly
those rows and the three columns used.

Full source (see DATA.md): regenerate ``immune_system_drugs.csv`` via
``src/data_build/build_immune_drug_dataset.py`` (Open Targets GraphQL API),
place it under data/raw/, then run this script.

Usage:
    python scripts/subset_approved_drugs.py \
        --input data/raw/immune_system_drugs.csv \
        --output data/reference/approved_target_drugs.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
KEEP = ["Target (Gene)", "Drug Name", "maxClinicalStage"]

DEFAULT_INPUT = REPO / "data" / "raw" / "immune_system_drugs.csv"
DEFAULT_OUTPUT = REPO / "data" / "reference" / "approved_target_drugs.csv"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    if not args.input.exists():
        raise SystemExit(
            f"raw drug pull not found: {args.input}\n"
            "Regenerate it via src/data_build/build_immune_drug_dataset.py (see DATA.md)."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    n = 0
    with open(args.input, newline="") as fh_in, open(args.output, "w", newline="") as fh_out:
        reader = csv.DictReader(fh_in)
        writer = csv.DictWriter(fh_out, fieldnames=KEEP)
        writer.writeheader()
        for row in reader:
            if row.get("maxClinicalStage") != "APPROVAL":
                continue
            key = (row["Target (Gene)"], row["Drug Name"])
            if key in seen:
                continue
            seen.add(key)
            writer.writerow({c: row[c] for c in KEEP})
            n += 1

    size_kb = args.output.stat().st_size / 1e3
    print(f"wrote {args.output} ({n} approved target-drug rows, {size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
