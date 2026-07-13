"""Combine the Spark `part-*.parquet` shards in this folder into a single CSV.

Nested columns (chemblIds, targets, references) are serialized as JSON strings
so the output is flat, lossless, and reloadable via json.loads.

Uses pyarrow + the stdlib csv module only (no pandas/numpy dependency).
"""

import csv
import json
from pathlib import Path

import pyarrow.parquet as pq

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "raw"
OUTPUT = DATA_DIR / "part_combined.csv"
NESTED_COLS = {"chemblIds", "targets", "references"}


def main():
    parts = sorted(DATA_DIR.glob("part-*.parquet"))
    if not parts:
        raise SystemExit(f"No part-*.parquet files found in {DATA_DIR}")

    columns = None
    total = 0
    with OUTPUT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        for part in parts:
            table = pq.ParquetFile(part).read()
            if columns is None:
                columns = table.column_names
                writer.writerow(columns)
            elif table.column_names != columns:
                raise SystemExit(
                    f"Schema mismatch in {part.name}: {table.column_names} != {columns}"
                )

            for row in table.to_pylist():
                writer.writerow(
                    json.dumps(row[col]) if col in NESTED_COLS else row[col]
                    for col in columns
                )
                total += 1

    print(f"Combined {len(parts)} part file(s) -> {OUTPUT}")
    print(f"Rows: {total}  Columns: {columns}")


if __name__ == "__main__":
    main()
