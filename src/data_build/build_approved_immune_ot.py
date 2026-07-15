"""Build the approved-immune-drug tables from Open Targets alone (no FDA layer).

Produces two CSVs, one row-level and one gene-level:

  approved_immune_drugs_evidence_ot.csv   -- one row per (drug x immune-indication)
  approved_immune_drugs_by_gene_ot.csv    -- gene rollup

Pipeline (Open Targets Platform GraphQL only):
    0. immune disease set  = descendants of MONDO_0005046 (immune system disorder)
    1. candidate drugs      = drugs hitting immune-associated targets that have an
                              APPROVAL-stage indication inside the immune set
    2. per-drug detail      = per-INDICATION maxClinicalStage, mechanism of action
                              (action type + text), gene target(s), parentMolecule
    3. salt-collapse        = re-map every salt/child to its parentMolecule so the
                              same active moiety is not counted twice
    4. rows                 = keep (drug x indication) where the indication is in the
                              immune set AND its per-indication stage == APPROVAL
    5. gene rollup

What "approved" means here (important):
    A row's ot_approved_immune == 1 means the (drug x indication) pair reached the
    APPROVAL stage FOR THAT SPECIFIC IMMUNE INDICATION, granted by at least one
    regulator SOMEWHERE in the world. Open Targets' clinical stage derives from
    ChEMBL max_phase and is regulator-agnostic: it does NOT distinguish FDA vs EMA
    vs PMDA, does not tell you the drug is still marketed, and carries no per-
    indication agency field (verified against the live schema). If you need
    US-specific ("FDA-approved for this indication") confirmation, that requires an
    external source such as openFDA SPL labels -- deliberately omitted from this
    OT-only build.

Design decisions worth knowing:
  * Unit is (drug x immune-indication), never per-drug: dual-use biologics
    (e.g. rituximab: RA/GPA/pemphigus AND lymphoma) keep every immune row and are
    never disqualified by their non-immune approvals.
  * OT `maximumClinicalStage` is the DRUG-level max across ALL indications; we use
    the PER-INDICATION stage instead, so "approved somewhere for something" is not
    mistaken for "approved for this immune disease".
  * Salt-collapse uses OT's own parentMolecule pointer (salt/parent share the active
    moiety, target, and MoA -- verified) so the same molecule is not counted twice.

Pure stdlib csv + requests (no pandas), matching build_immune_drug_dataset.py.

Usage:
    python build_approved_immune_ot.py \
        --evidence-out approved_immune_drugs_evidence_ot.csv \
        --gene-out     approved_immune_drugs_by_gene_ot.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import requests

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
PARENT_ID = "MONDO_0005046"  # immune system disorder
PAGE = 100
BATCH = 20  # drugs per aliased GraphQL query

EVIDENCE_COLS = [
    "drug", "chembl_id", "drug_type", "gene_target", "action_type", "moa",
    "immune_indication", "disease_id", "drug_stage_for_immune_indication",
    "ot_approved_immune",
]
GENE_COLS = [
    "gene_target", "n_approved_drugs", "action_types", "approved_drugs",
    "immune_indications",
]


# --------------------------------------------------------------------------- #
# transport                                                                   #
# --------------------------------------------------------------------------- #
def ot_gql(session, query, variables=None, retries=5):
    """POST a GraphQL query to Open Targets with retry/backoff."""
    for attempt in range(retries):
        try:
            resp = session.post(OT_URL, json={"query": query, "variables": variables or {}},
                                timeout=120)
            if resp.status_code == 200:
                payload = resp.json()
                if payload.get("errors"):
                    raise RuntimeError(payload["errors"])
                return payload["data"]
            last = f"HTTP {resp.status_code}: {resp.text[:200]}"
        except Exception as exc:  # network / JSON / API error
            last = str(exc)
        wait = 2 ** attempt
        print(f"  OT retry {attempt + 1}/{retries} ({last[:100]}); sleeping {wait}s",
              file=sys.stderr)
        time.sleep(wait)
    raise RuntimeError(f"OT query failed after {retries} attempts: {last}")


# --------------------------------------------------------------------------- #
# step 0 - immune disease id set                                              #
# --------------------------------------------------------------------------- #
def immune_disease_ids(session):
    q = "query($e:String!){disease(efoId:$e){descendants}}"
    d = ot_gql(session, q, {"e": PARENT_ID})["disease"]
    return set(d["descendants"]) | {PARENT_ID}


# --------------------------------------------------------------------------- #
# step 1 - candidate drugs (approval-stage, immune indication)                #
# --------------------------------------------------------------------------- #
CAND_Q = """
query($e:String!,$i:Int!,$s:Int!){
  disease(efoId:$e){
    associatedTargets(enableIndirect:true, page:{index:$i,size:$s}){
      count
      rows{ target{ approvedSymbol id
        drugAndClinicalCandidates{ rows{ maxClinicalStage
          drug{ id name }
          diseases{ disease{ id } } } } } }
    }
  }
}"""


def candidate_drug_ids(session, immune_ids):
    drug_ids = {}
    index, total = 0, None
    while True:
        at = ot_gql(session, CAND_Q, {"e": PARENT_ID, "i": index, "s": PAGE})[
            "disease"]["associatedTargets"]
        if total is None:
            total = at["count"]
            print(f"  associated targets: {total}")
        rows = at["rows"]
        if not rows:
            break
        for row in rows:
            for cand in row["target"]["drugAndClinicalCandidates"]["rows"]:
                drug = cand.get("drug") or {}
                did = drug.get("id")
                if not did or cand.get("maxClinicalStage") != "APPROVAL":
                    continue
                dis = {d["disease"]["id"] for d in (cand.get("diseases") or [])
                       if d.get("disease")}
                if dis & immune_ids:
                    drug_ids[did] = drug.get("name")
        index += 1
        print(f"    page {index} | candidates so far: {len(drug_ids)}")
        if index * PAGE >= total:
            break
    return drug_ids


# --------------------------------------------------------------------------- #
# step 2 - per-drug detail (aliased batches)                                  #
# --------------------------------------------------------------------------- #
DETAIL_FRAG = """
  d%d: drug(chemblId:"%s"){
    id name drugType maximumClinicalStage
    parentMolecule{ id name }
    mechanismsOfAction{ rows{ mechanismOfAction actionType targets{ approvedSymbol id } } }
    indications{ rows{ maxClinicalStage disease{ id name } } }
  }"""


def fetch_detail(session, chembl_ids):
    detail = {}
    ids = list(chembl_ids)
    for i in range(0, len(ids), BATCH):
        batch = ids[i:i + BATCH]
        q = "{" + "".join(DETAIL_FRAG % (j, cid) for j, cid in enumerate(batch)) + "}"
        data = ot_gql(session, q)
        for j, cid in enumerate(batch):
            detail[cid] = data.get(f"d{j}")
        if (i // BATCH) % 10 == 0:
            print(f"    detail {min(i + BATCH, len(ids))}/{len(ids)}")
    return detail


# --------------------------------------------------------------------------- #
# step 3 - salt-collapse via parentMolecule                                   #
# --------------------------------------------------------------------------- #
def collapse_to_parents(session, detail):
    """Map each chembl id to its canonical parent; backfill missing parent detail."""
    parent_of, name_of = {}, {}
    for cid, d in detail.items():
        if not d:
            parent_of[cid] = cid
            continue
        pm = d.get("parentMolecule")
        parent_of[cid] = pm["id"] if pm else cid
        name_of[cid] = d.get("name")
        if pm:
            name_of[pm["id"]] = pm["name"]
    missing = [cid for cid in set(parent_of.values())
               if cid not in detail or not detail.get(cid)]
    if missing:
        print(f"  backfilling detail for {len(missing)} parent molecules")
        detail.update(fetch_detail(session, missing))
    return parent_of, name_of, detail


# --------------------------------------------------------------------------- #
# step 4 - (drug x immune-indication) rows on canonical parents               #
# --------------------------------------------------------------------------- #
def build_rows(detail, parent_of, immune_ids):
    rows, seen = [], set()
    for cid in sorted(set(parent_of.values())):
        d = detail.get(cid)
        if not d:
            continue
        moa = d.get("mechanismsOfAction", {}).get("rows") or []
        genes = sorted({t["approvedSymbol"] for m in moa for t in (m.get("targets") or [])})
        acts = sorted({m.get("actionType") for m in moa if m.get("actionType")})
        moatext = sorted({m.get("mechanismOfAction") for m in moa if m.get("mechanismOfAction")})
        for ind in (d.get("indications", {}).get("rows") or []):
            dis = ind.get("disease") or {}
            if dis.get("id") in immune_ids and ind.get("maxClinicalStage") == "APPROVAL":
                key = (cid, dis.get("id"))
                if key in seen:
                    continue
                seen.add(key)
                rows.append({
                    "chembl_id": cid, "drug": d.get("name"), "drug_type": d.get("drugType"),
                    "gene_target": ";".join(genes), "action_type": ";".join(acts),
                    "moa": ";".join(moatext),
                    "immune_indication": dis.get("name"), "disease_id": dis.get("id"),
                    "drug_stage_for_immune_indication": ind.get("maxClinicalStage"),
                    "ot_approved_immune": 1,
                })
    return rows


# --------------------------------------------------------------------------- #
# step 5 - gene rollup                                                        #
# --------------------------------------------------------------------------- #
def gene_rollup(rows):
    gene = defaultdict(lambda: {"drugs": set(), "acts": set(), "inds": set()})
    for r in rows:
        for g in r["gene_target"].split(";"):
            if not g:
                continue
            gene[g]["drugs"].add(r["drug"])
            # r["action_type"] is a ";"-joined string; split so acts holds single tokens
            for a in r["action_type"].split(";"):
                if a:
                    gene[g]["acts"].add(a)
            gene[g]["inds"].add(r["immune_indication"])
    return [{
        "gene_target": g,
        "n_approved_drugs": len(v["drugs"]),
        "action_types": ";".join(sorted(set(v["acts"]))),
        "approved_drugs": ";".join(sorted(set(v["drugs"]))[:15]),
        "immune_indications": ";".join(sorted(set(v["inds"]))[:12]),
    } for g, v in sorted(gene.items(), key=lambda kv: -len(kv[1]["drugs"]))]


def write_csv(path, cols, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence-out", default="approved_immune_drugs_evidence_ot.csv")
    ap.add_argument("--gene-out", default="approved_immune_drugs_by_gene_ot.csv")
    args = ap.parse_args()

    session = requests.Session()

    print("0. immune disease id set ...")
    immune_ids = immune_disease_ids(session)
    print(f"   immune ids: {len(immune_ids)}")

    print("1. candidate drugs ...")
    drug_ids = candidate_drug_ids(session, immune_ids)
    print(f"   candidate drugs: {len(drug_ids)}")

    print("2. per-drug detail ...")
    detail = fetch_detail(session, drug_ids)

    print("3. salt-collapse via parentMolecule ...")
    parent_of, name_of, detail = collapse_to_parents(session, detail)
    n_salt = sum(1 for k, v in parent_of.items() if k != v)
    print(f"   salts re-mapped: {n_salt} | canonical molecules: {len(set(parent_of.values()))}")

    print("4. (drug x immune-indication) rows ...")
    rows = build_rows(detail, parent_of, immune_ids)
    rows.sort(key=lambda r: (r["drug"], r["immune_indication"]))
    print(f"   rows: {len(rows)} | drugs: {len({r['chembl_id'] for r in rows})}")

    print("5. gene rollup ...")
    grows = gene_rollup(rows)

    write_csv(args.evidence_out, EVIDENCE_COLS, rows)
    write_csv(args.gene_out, GENE_COLS, grows)
    print(f"\nwrote {args.evidence_out} ({len(rows)} rows) and "
          f"{args.gene_out} ({len(grows)} genes)")
    print(f"  drug_type mix: {dict(Counter(r['drug_type'] for r in rows))}")


if __name__ == "__main__":
    main()
