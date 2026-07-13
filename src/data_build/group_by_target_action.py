"""Collapse the approved-immune dataset to one row per (Target, Action Type).

Rows sharing the same Target (Gene) / Action Type are merged; their distinct Actions
(mechanism of action), Drug Names, and Specific Diseases are each combined into
"; "-joined lists (semicolon, not comma, because many disease names contain internal
commas). Adds num_drugs / num_diseases counts.
"""

import csv
from collections import OrderedDict
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
SOURCE = _REPO / "data" / "raw" / "immune_system_drugs_immune_approved.csv"
OUTPUT = _REPO / "data" / "derived" / "immune_system_drugs_grouped.csv"
SEP = "; "

KEY_COLS = ["Target (Gene)", "Action Type"]
# Columns collapsed into distinct "; "-joined lists (Action is no longer a group key).
LIST_COLS = ["Action", "Drug Name", "Specific Disease"]
CONST_COLS = ["Broad Category", "maxClinicalStage", "immune_indication"]
FIELDNAMES = KEY_COLS + LIST_COLS + ["num_drugs", "num_diseases"] + CONST_COLS


def main():
    # group key -> {"lists": {col: OrderedDict-as-set}, "const": {col: value}}
    groups = OrderedDict()
    total = 0
    with SOURCE.open(newline="", encoding="utf-8") as fin:
        for row in csv.DictReader(fin):
            total += 1
            key = tuple(row[c] for c in KEY_COLS)
            grp = groups.get(key)
            if grp is None:
                grp = {"lists": {c: OrderedDict() for c in LIST_COLS},
                       "const": {c: row[c] for c in CONST_COLS}}
                groups[key] = grp
            for c in LIST_COLS:
                grp["lists"][c].setdefault(row[c], None)
            # defensively confirm constant columns really are constant per group
            for c in CONST_COLS:
                if row[c] != grp["const"][c]:
                    raise ValueError(
                        f"Column {c!r} not constant for group {key}: "
                        f"{row[c]!r} vs {grp['const'][c]!r}"
                    )

    with OUTPUT.open("w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=FIELDNAMES)
        writer.writeheader()
        for key, grp in groups.items():
            out = dict(zip(KEY_COLS, key))
            for c in LIST_COLS:
                out[c] = SEP.join(sorted(v for v in grp["lists"][c] if v))
            out["num_drugs"] = len([v for v in grp["lists"]["Drug Name"] if v])
            out["num_diseases"] = len([v for v in grp["lists"]["Specific Disease"] if v])
            out.update(grp["const"])
            writer.writerow(out)

    print(f"Read {total} rows from {SOURCE.name}")
    print(f"Wrote {len(groups)} grouped rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
