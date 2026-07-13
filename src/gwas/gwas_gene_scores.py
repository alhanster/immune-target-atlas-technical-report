"""
Query Open Targets for the GWAS (gwas_credible_sets) gene-level
association score for every target associated with a disease.

Datasource `gwas_credible_sets` is the GWAS component of the
`genetic_association` datatype. The `datasources` filter restricts the
result to targets that carry GWAS evidence; the per-target GWAS score is
read from `datasourceScores` (id == "gwas_credible_sets") and the table
is sorted on it. Requires only `requests`.

Outputs (written to data/derived/ when run as a script):
  - gwas_gene_scores_<efo>.csv: rank, gene, ensembl_id, gwas_score for
    the GWAS-evidenced genes only (harmonized to current HGNC).
  - gwas_gene_scores.csv: the full gnomad gene universe (gene column only,
    harmonized) left-joined with gwas_score, with 0 for any gnomad gene lacking
    GWAS evidence. Consumed by src/data_build/creating_full_gene_list.py.
"""
import requests

API = "https://api.platform.opentargets.org/api/v4/graphql"

QUERY = """
query($efo: String!, $size: Int!, $index: Int!) {
  disease(efoId: $efo) {
    id
    name
    associatedTargets(
      page: { size: $size, index: $index }
      datasources: [
        { id: "gwas_credible_sets", weight: 1.0, required: true, propagate: true }
      ]
    ) {
      count
      rows {
        target { id approvedSymbol }
        datasourceScores { id score }
      }
    }
  }
}
"""


def gwas_gene_scores(efo_id, page_size=200, session=None):
    """Yield (rank-agnostic) dicts of GWAS gene-level scores for a disease.

    Pages through all associated targets and returns one record per gene:
    {symbol, ensembl_id, gwas_score}.
    """
    s = session or requests.Session()
    index, total, out = 0, None, []
    while True:
        variables = {"efo": efo_id, "size": page_size, "index": index}
        resp = s.post(API, json={"query": QUERY, "variables": variables}, timeout=60)
        resp.raise_for_status()
        payload = resp.json()
        if payload.get("errors"):
            raise RuntimeError(payload["errors"])
        disease = payload["data"]["disease"]
        if disease is None:
            raise ValueError(f"Unknown disease id: {efo_id}")
        at = disease["associatedTargets"]
        total = at["count"]
        rows = at["rows"]
        if not rows:
            break
        for row in rows:
            gwas = next(
                (d["score"] for d in row["datasourceScores"]
                 if d["id"] == "gwas_credible_sets"),
                None,
            )
            out.append({
                "symbol": row["target"]["approvedSymbol"],
                "ensembl_id": row["target"]["id"],
                "gwas_score": gwas,
            })
        index += 1
        if len(out) >= total:
            break
    out.sort(key=lambda r: (r["gwas_score"] is None, -(r["gwas_score"] or 0)))
    return disease["name"], total, out


if __name__ == "__main__":
    import sys
    from pathlib import Path

    import pandas as pd

    # This script lives in src/gwas/; shared gene-name utilities live in
    # src/shared/ and derived outputs go to data/derived/.
    HERE = Path(__file__).resolve().parent        # src/gwas
    REPO = HERE.parents[1]                          # repo root
    DERIVED = REPO / "data" / "derived"
    sys.path.insert(0, str(REPO / "src" / "shared"))
    from gene_name_utils import load_hgnc, harmonize, check_format

    efo = "MONDO_0005046"
    name, total, recs = gwas_gene_scores(efo)
    print(f"{efo}  ({name}): {total} targets with GWAS evidence")

    df = pd.DataFrame(recs).rename(columns={"symbol": "gene"})

    # Harmonize gene symbols to current HGNC + formatting check before saving.
    hgnc = load_hgnc()                                   # loaded once, reused below
    df, _ = harmonize(df, hgnc, gene_col="gene")
    check_format(df["gene"])

    # Records are already sorted by gwas_score (desc) and harmonize preserves
    # order; assign a contiguous rank after harmonization.
    df.insert(0, "rank", range(1, len(df) + 1))
    df = df[["rank", "gene", "ensembl_id", "gwas_score"]]

    # Full table (GWAS-evidenced genes only), ranked.
    out_csv = DERIVED / f"gwas_gene_scores_{efo}.csv"
    df.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"wrote {len(df)} rows -> {out_csv}")
    print(df.head(15).to_string(index=False))

    # Slim [gene, gwas_score] table over the FULL gnomad gene universe: take only
    # gnomad's gene column, harmonize it to current HGNC (so genes gnomad lists
    # under old symbols still match the harmonized GWAS scores), left-join the
    # score, and write 0 for any gnomad gene with no GWAS evidence.
    gnomad = pd.read_csv(REPO / "data" / "reference" / "gnomad_constraint_subset.tsv",
                         sep="\t", dtype=str, keep_default_na=False, usecols=["gene"])
    gnomad = gnomad[(gnomad["gene"].str.strip() != "") & (gnomad["gene"] != "NA")]
    gnomad = gnomad.drop_duplicates("gene")
    gnomad, _ = harmonize(gnomad, hgnc, gene_col="gene")
    slim = gnomad.merge(df[["gene", "gwas_score"]], on="gene", how="left")
    slim["gwas_score"] = slim["gwas_score"].fillna(0)
    slim = slim.sort_values("gwas_score", ascending=False)

    # Written to data/derived/ (consumed by creating_full_gene_list.py).
    for dest in (DERIVED / "gwas_gene_scores.csv",):
        slim.to_csv(dest, index=False, float_format="%.6f")
        print(f"wrote {len(slim)} rows ({int((slim['gwas_score'] > 0).sum())} with a "
              f"GWAS score) -> {dest}")
