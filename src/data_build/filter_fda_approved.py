"""Subset of immune_system_drugs.csv restricted to FDA-approved drugs.

"FDA-approved" maps to Open Targets `maxClinicalStage == "APPROVAL"` (the highest
clinical stage; there is no integer `phase` field in the current schema).
"""

import csv
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "raw"
SOURCE = DATA_DIR / "immune_system_drugs.csv"
OUTPUT = DATA_DIR / "immune_system_drugs_fda_approved.csv"
APPROVED_STAGE = "APPROVAL"


def main():
    with SOURCE.open(newline="", encoding="utf-8") as fin, \
            OUTPUT.open("w", newline="", encoding="utf-8") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
        writer.writeheader()
        kept = total = 0
        for row in reader:
            total += 1
            if row["maxClinicalStage"] == APPROVED_STAGE:
                writer.writerow(row)
                kept += 1

    print(f"Read {total} rows from {SOURCE.name}")
    print(f"Kept {kept} FDA-approved (maxClinicalStage == {APPROVED_STAGE}) -> {OUTPUT}")


if __name__ == "__main__":
    main()
