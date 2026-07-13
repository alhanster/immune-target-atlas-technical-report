"""Left-join the PU score onto the full gene list.

Produces `scored_full_gene_list`: every gene in `data/derived/full_gene_list.tsv`
with its `pu_score` (from the pipeline's `outputs/nominations.tsv`) inserted right
after `gene_id`. All original columns and rows are preserved (left join).

    python src/model/build_scored_full_gene_list.py

Reads  : data/derived/full_gene_list.tsv, outputs/nominations.tsv
Writes : outputs/scored_full_gene_list.tsv
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

# build_scored_full_gene_list.py lives at <repo>/src/model/ -> parents[2] is repo root.
ROOT = Path(__file__).resolve().parents[2]
GENE_LIST_TSV = ROOT / "data" / "derived" / "full_gene_list.tsv"
NOMINATIONS_TSV = ROOT / "outputs" / "nominations.tsv"
OUT_TSV = ROOT / "outputs" / "scored_full_gene_list.tsv"


def main() -> None:
    full_gene_list = pd.read_csv(GENE_LIST_TSV, sep="\t")
    if not NOMINATIONS_TSV.exists():
        raise SystemExit(
            f"{NOMINATIONS_TSV} not found — run `python -m model.run` first "
            "to generate the PU scores."
        )
    nom = pd.read_csv(NOMINATIONS_TSV, sep="\t")

    # Left join on the stable Ensembl id: keep every gene-list row and column,
    # attach pu_score.
    scored_full_gene_list = full_gene_list.merge(
        nom[["gene_id", "pu_score"]], on="gene_id", how="left"
    )

    # Move pu_score to sit immediately after gene_id.
    cols = list(scored_full_gene_list.columns)
    cols.remove("pu_score")
    cols.insert(cols.index("gene_id") + 1, "pu_score")
    scored_full_gene_list = scored_full_gene_list[cols]

    # Guardrails: no rows dropped, every gene scored.
    assert len(scored_full_gene_list) == len(full_gene_list), "row count changed"
    missing = int(scored_full_gene_list["pu_score"].isna().sum())
    assert missing == 0, (
        f"{missing} genes have no pu_score — nominations.tsv may be stale; "
        "re-run `python -m model.run`."
    )

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    scored_full_gene_list.to_csv(
        OUT_TSV, sep="\t", index=False, float_format="%.6f"
    )

    print(f"scored_full_gene_list.shape = {scored_full_gene_list.shape}")
    print(f"columns: {list(scored_full_gene_list.columns)}")
    print(scored_full_gene_list.head().to_string(index=False))
    print(f"\nwrote {OUT_TSV}")


if __name__ == "__main__":
    main()
