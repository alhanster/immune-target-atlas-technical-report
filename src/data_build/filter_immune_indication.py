"""Subset of immune_system_drugs.csv restricted to immune-system indications.

Keeps rows where immune_indication == 1, i.e. the drug's indication is
MONDO_0005046 ("immune system disorder") or one of its descendants.
Clinical stage is NOT filtered here (all stages kept).
"""

import csv
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "raw"
SOURCE = DATA_DIR / "immune_system_drugs.csv"
OUTPUT = DATA_DIR / "immune_system_drugs_immune_only.csv"


def main():
    drugs = set()
    with SOURCE.open(newline="", encoding="utf-8") as fin, \
            OUTPUT.open("w", newline="", encoding="utf-8") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
        writer.writeheader()
        kept = total = 0
        for row in reader:
            total += 1
            if row["immune_indication"] == "1":
                writer.writerow(row)
                drugs.add(row["Drug Name"])
                kept += 1

    print(f"Read {total} rows from {SOURCE.name}")
    print(f"Kept {kept} immune-indication rows ({len(drugs)} distinct drugs) -> {OUTPUT}")


if __name__ == "__main__":
    main()
