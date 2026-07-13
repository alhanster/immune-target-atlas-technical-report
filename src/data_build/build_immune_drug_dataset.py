"""Build a drug dataset for immune-system diseases from the Open Targets GraphQL API.

This is a corrected rewrite of an older script that targeted a now-removed schema.
Current schema facts (verified against the live API, 2026-07):
  - Disease IDs are MONDO, not EFO. Immune system disorder = MONDO_0005046.
  - `associatedTargets` uses `enableIndirect` (not `enableDirect`).
  - `Target.knownDrugs` is gone; use `Target.drugAndClinicalCandidates`, whose rows carry
    `maxClinicalStage` (enum string, e.g. APPROVAL / PHASE_3), the `drug` (with
    `mechanismsOfAction`), and the drug's indications in `diseases[]`.

Output: Data/immune_system_drugs.csv (one row per target x drug-candidate x indication).
Pure stdlib csv + requests (no pandas: Anaconda's pandas/numpy is binary-broken here).
"""

import csv
import sys
import time
from pathlib import Path

import requests

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
PARENT_ID = "MONDO_0005046"  # immune system disorder
BROAD_CATEGORY = "immune system disorder"
PAGE = 100
OUTPUT = Path(__file__).resolve().parents[2] / "data" / "raw" / "immune_system_drugs.csv"

FIELDNAMES = [
    "Target (Gene)",
    "Drug Name",
    "Action",
    "Action Type",
    "Specific Disease",
    "Specific Disease Id",
    "Broad Category",
    "maxClinicalStage",
    "immune_indication",
]

DESCENDANTS_QUERY = """
query($efoId: String!) {
  disease(efoId: $efoId) { id descendants }
}
"""

TARGETS_QUERY = """
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
              drug {
                name
                mechanismsOfAction { rows { mechanismOfAction actionType } }
              }
              diseases { diseaseFromSource disease { id name } }
            }
          }
        }
      }
    }
  }
}
"""


def gql(session, query, variables, retries=5):
    """POST a GraphQL query with retry/backoff on transport or API errors."""
    for attempt in range(retries):
        try:
            resp = session.post(
                API_URL, json={"query": query, "variables": variables}, timeout=120
            )
            if resp.status_code == 200:
                payload = resp.json()
                if payload.get("errors"):
                    raise RuntimeError(payload["errors"])
                return payload["data"]
            last = f"HTTP {resp.status_code}: {resp.text[:200]}"
        except Exception as exc:  # network, JSON, or API errors
            last = str(exc)
        wait = 2 ** attempt
        print(f"  retry {attempt + 1}/{retries} after error ({last[:120]}); "
              f"sleeping {wait}s", file=sys.stderr)
        time.sleep(wait)
    raise RuntimeError(f"Query failed after {retries} attempts: {last}")


def join_distinct(values):
    """Join distinct, order-preserving, non-empty strings with '; '."""
    seen = []
    for v in values:
        if v and v not in seen:
            seen.append(v)
    return "; ".join(seen)


def build_rows(target, immune_ids):
    """Yield flattened dataset rows for one target node."""
    symbol = target["approvedSymbol"]
    for cand in target["drugAndClinicalCandidates"]["rows"]:
        drug = cand.get("drug") or {}
        moa_rows = (drug.get("mechanismsOfAction") or {}).get("rows") or []
        action = join_distinct(m.get("mechanismOfAction") for m in moa_rows)
        action_type = join_distinct(m.get("actionType") for m in moa_rows)
        base = {
            "Target (Gene)": symbol,
            "Drug Name": drug.get("name"),
            "Action": action,
            "Action Type": action_type,
            "Broad Category": BROAD_CATEGORY,
            "maxClinicalStage": cand.get("maxClinicalStage"),
        }
        indications = cand.get("diseases") or []
        if not indications:
            yield {**base, "Specific Disease": "", "Specific Disease Id": "",
                   "immune_indication": 0}
            continue
        for ind in indications:
            dis = ind.get("disease") or {}
            dis_id = dis.get("id") or ""
            name = dis.get("name") or ind.get("diseaseFromSource") or ""
            yield {**base,
                   "Specific Disease": name,
                   "Specific Disease Id": dis_id,
                   "immune_indication": 1 if dis_id in immune_ids else 0}


def main():
    session = requests.Session()

    # 1. Descendant set (for the immune_indication flag).
    print(f"Fetching descendants of {PARENT_ID} ...")
    disease = gql(session, DESCENDANTS_QUERY, {"efoId": PARENT_ID})["disease"]
    immune_ids = set(disease["descendants"]) | {PARENT_ID}
    print(f"  immune-disorder id set size: {len(immune_ids)}")

    # 2/3. Paginate associated targets, flatten to rows, de-dup.
    seen = set()
    rows = []
    index = 0
    total = None
    while True:
        data = gql(session, TARGETS_QUERY,
                   {"efoId": PARENT_ID, "index": index, "size": PAGE})
        at = data["disease"]["associatedTargets"]
        if total is None:
            total = at["count"]
            print(f"Associated targets: {total} (page size {PAGE}, "
                  f"~{-(-total // PAGE)} pages)")
        page_rows = at["rows"]
        if not page_rows:
            break
        for row in page_rows:
            for out in build_rows(row["target"], immune_ids):
                key = tuple(out[f] for f in FIELDNAMES)
                if key not in seen:
                    seen.add(key)
                    rows.append(out)
        index += 1
        print(f"  page {index} done | targets scanned: {index * PAGE} | "
              f"rows so far: {len(rows)}")
        if index * PAGE >= total:
            break

    # 4. Write CSV.
    with OUTPUT.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    n_immune = sum(r["immune_indication"] for r in rows)
    print(f"\nWrote {len(rows)} rows -> {OUTPUT}")
    print(f"  immune_indication=1: {n_immune} | =0: {len(rows) - n_immune}")


if __name__ == "__main__":
    main()
