# How `immune_system_drugs_grouped.csv` is produced — query steps

This documents the full, reproducible pipeline from the Open Targets GraphQL API to the
final grouped file. Each step names its script, input, output, and row count.

Source API: `https://api.platform.opentargets.org/api/v4/graphql` (Open Targets Platform, v4).
All scripts use `requests` + the Python stdlib `csv` module (no pandas). Run them from the
project root, e.g. `python3 Data/<script>.py`.

Pipeline overview:

```
Open Targets GraphQL  ──①──▶ immune_system_drugs.csv           (323,855 rows)
                                     │
                                     ├──② immune_indication==1 ─▶ immune_system_drugs_immune_only.csv     (21,735)
                                     │
        immune_only ────────────────┴──③ maxClinicalStage==APPROVAL ─▶ immune_system_drugs_immune_approved.csv (18,202)
                                                                                 │
                                                                                 └──④ group by Target + Action Type
                                                                                        └─▶ immune_system_drugs_grouped.csv (983)
```

---

## Step 1 — Query Open Targets → `immune_system_drugs.csv`

Script: `Data/build_immune_drug_dataset.py`. Output: 323,855 rows.

**1a. Get the immune-disease id set.** One query resolves every disease under
"immune system disorder" (`MONDO_0005046`). `descendants` is a native Open Targets field
returning a flat list of MONDO id strings (the ontology sub-tree):

```graphql
query($efoId: String!) {
  disease(efoId: $efoId) { id descendants }
}
```
Variables: `{ "efoId": "MONDO_0005046" }`.
In code the returned list is turned into a set and the parent id is added:
`immune_ids = set(descendants) | {"MONDO_0005046"}` → **1,189 ids**. Used only for the
`immune_indication` flag in step 1c.

**1b. Page through the disease's associated targets and their drugs.** The parent disease has
19,067 indirectly-associated targets; each is fetched with its drug candidates nested inside a
paginated query (`page.size = 100`, ~191 pages):

```graphql
query($efoId: String!, $index: Int!, $size: Int!) {
  disease(efoId: $efoId) {
    associatedTargets(enableIndirect: true, page: {index: $index, size: $size}) {
      count
      rows {
        target {
          approvedSymbol
          drugAndClinicalCandidates {
            rows {
              maxClinicalStage
              drug { name mechanismsOfAction { rows { mechanismOfAction actionType } } }
              diseases { diseaseFromSource disease { id name } }
            }
          }
        }
      }
    }
  }
}
```

**1c. Flatten to one row per (target × drug-candidate × indication).** Columns written:

| Column | From |
|---|---|
| `Target (Gene)` | `target.approvedSymbol` |
| `Drug Name` | `drug.name` |
| `Action` | distinct `drug.mechanismsOfAction.rows[].mechanismOfAction`, joined `"; "` |
| `Action Type` | distinct `drug.mechanismsOfAction.rows[].actionType`, joined `"; "` |
| `Specific Disease` | `diseases[].disease.name` (fallback `diseaseFromSource`) |
| `Specific Disease Id` | `diseases[].disease.id` |
| `Broad Category` | constant `"immune system disorder"` |
| `maxClinicalStage` | `drugAndClinicalCandidates.rows[].maxClinicalStage` (enum, e.g. `APPROVAL`, `PHASE_3`) |
| `immune_indication` | `1` if `Specific Disease Id ∈ immune_ids` else `0` (derived — see step 1a) |

Exact-duplicate rows are dropped.

> **Schema note (why this differs from the original script):** the old query used
> `EFO_0000540` (now migrated to MONDO), `enableDirect` (now `enableIndirect`), and
> `target.knownDrugs { … phase }` (removed; replaced by `drugAndClinicalCandidates` with
> `maxClinicalStage` — there is no integer `phase`).

---

## Step 2 — Keep immune indications → `immune_system_drugs_immune_only.csv`

Script: `Data/filter_immune_indication.py`. Keeps rows where `immune_indication == "1"`.
323,855 → **21,735 rows** (1,611 distinct drugs). All clinical stages retained.

---

## Step 3 — Keep FDA-approved → `immune_system_drugs_immune_approved.csv`

Script: chained filter on step 2 for `maxClinicalStage == "APPROVAL"` (the current-schema
stand-in for "FDA-approved" / old Phase 4). 21,735 → **18,202 rows** (927 distinct drugs).
Every row now has `immune_indication == 1` **and** `maxClinicalStage == APPROVAL`.

---

## Step 4 — Group → `immune_system_drugs_grouped.csv`

Script: `Data/group_by_target_action.py`. **Group key: `Target (Gene)` + `Action Type`.**
18,202 → **983 rows**.

For each group, these columns are collapsed to distinct, alphabetically-sorted, `"; "`-joined
lists: `Action`, `Drug Name`, `Specific Disease`. Added: `num_drugs`, `num_diseases` (distinct
counts). Constant columns carried through as single values: `Broad Category`,
`maxClinicalStage` (`APPROVAL`), `immune_indication` (`1`). `Specific Disease Id` is dropped.

Separator is `"; "` (semicolon), **not** comma, because 423 disease names contain internal
commas (e.g. `"chronic myelogenous leukemia, BCR-ABL1 positive"`) — so a comma delimiter could
not be split back reliably.

Final columns:
`Target (Gene)`, `Action Type`, `Action`, `Drug Name`, `Specific Disease`, `num_drugs`,
`num_diseases`, `Broad Category`, `maxClinicalStage`, `immune_indication`.

> **Caveat:** `Action`, `Drug Name`, and `Specific Disease` are three *independent* lists within
> a group — the per-drug→per-disease pairing and per-drug mechanism are not preserved. Also,
> because a drug can hit multiple targets, a drug's non-grouped-target mechanisms may appear in
> the `Action` list. For row-level detail, use `immune_system_drugs_immune_approved.csv`.

---

## How to query / load the grouped file

```python
import csv

with open("Data/immune_system_drugs_grouped.csv", newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

# split a combined column back into a list
drugs = rows[0]["Drug Name"].split("; ")

# find every group targeting a gene
btk = [r for r in rows if r["Target (Gene)"] == "BTK"]

# groups with the most distinct drugs
top = sorted(rows, key=lambda r: int(r["num_drugs"]), reverse=True)[:10]
```

With pandas (if the environment's pandas is working):
```python
import pandas as pd
df = pd.read_csv("Data/immune_system_drugs_grouped.csv")
df["drug_list"] = df["Drug Name"].str.split("; ")
df.sort_values("num_drugs", ascending=False).head(10)
```
