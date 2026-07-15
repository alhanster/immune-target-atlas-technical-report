"""Data loading and gene-family grouping.

Loads the full gene table (all 18,692 rows kept — never dropped), the
approved-target list, and derives a per-gene CV family from HGNC gene groups
(primary) with a symbol-prefix fallback (approximate) for genes lacking a group.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import config


def load_gene_table(path: Path = config.GENE_LIST_TSV) -> pd.DataFrame:
    """Load full_gene_list.tsv. NA strings become real NaN; no rows dropped."""
    df = pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=True)
    df["gene"] = df["gene"].astype(str)
    # gwas_score is already 0-imputed upstream; guard against any stray NA.
    df["gwas_score"] = pd.to_numeric(df["gwas_score"], errors="coerce").fillna(0.0)
    df["IEI"] = pd.to_numeric(df["IEI"], errors="coerce").fillna(0).astype(int)
    return df


def load_approved_targets(path: Path = config.APPROVED_TXT) -> set[str]:
    """Return the set of approved drug-target gene symbols."""
    with open(path) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def load_family_groups(
    genes: list[str], path: Path = config.HGNC_GROUPS_TSV
) -> tuple[pd.Series, dict]:
    """Map each gene to a CV family id.

    Primary source: the FIRST HGNC ``gene_group_id`` for the gene. HGNC groups
    are the task's first-preference grouping. Connected-components over shared
    groups were rejected (a hub-group giant component swallows ~60% of genes);
    first-group-id gives a sane distribution (max family ~4% of genes).

    Fallback (flagged approximate): genes with no HGNC group get a coarse family
    from the 3-char gene-symbol prefix.

    Returns (family Series aligned to ``genes``, info dict).
    """
    fam_by_gene: dict[str, str] = {}
    if path.exists():
        hg = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        for sym, gid in zip(hg["symbol"], hg["gene_group_id"]):
            toks = [t for t in str(gid).split("|") if t.strip()]
            if toks:
                fam_by_gene[sym] = "HGNC:" + toks[0]

    families, n_hgnc, n_prefix = [], 0, 0
    for g in genes:
        if g in fam_by_gene:
            families.append(fam_by_gene[g]); n_hgnc += 1
        else:
            families.append("PFX:" + g[:3].upper()); n_prefix += 1

    fam = pd.Series(families, index=genes, name="family")
    info = {
        "grouping_source": "HGNC gene_group_id (first token); symbol-prefix fallback",
        "n_genes": len(genes),
        "n_families": int(fam.nunique()),
        "n_from_hgnc": n_hgnc,
        "n_from_prefix_fallback": n_prefix,
        "prefix_fallback_is_approximate": True,
        "hgnc_reference_present": path.exists(),
        "largest_family_frac": float(fam.value_counts().iloc[0] / len(genes)),
    }
    return fam, info
