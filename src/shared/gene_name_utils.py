#!/usr/bin/env python3
"""Shared gene-name utilities: formatting checks + HGNC harmonization.

One place for the gene-name logic reused across the scripts in this directory,
so formatting rules and HGNC symbol handling stay consistent.

The HGNC reference is fetched from the HGNC complete set and cached, so rename
maps are derived dynamically (no hardcoded, drifting dicts). Point `load_hgnc`
at a local TSV to stay offline.

Functions
    load_hgnc(path=None, refresh=False) -> dict(approved, prev2app, alias)
    check_format(genes, verbose=True)   -> dict report (no network)
    classify_hgnc(genes, hgnc, verbose=True) -> pd.Series of statuses
    harmonize(df, hgnc, gene_col="gene", verbose=True) -> (df, info)

CLI
    python gene_name_utils.py <csv> [gene_col]   # format + HGNC report
"""

from __future__ import annotations

import re
import sys
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd

HGNC_URL = ("https://storage.googleapis.com/public-download-files/hgnc/tsv/"
            "tsv/hgnc_complete_set.txt")
HGNC_CACHE = Path(tempfile.gettempdir()) / "hgnc_complete_set.txt"

# Approved HGNC symbols are uppercase alnum plus a little punctuation; the body
# may be lowercase for the `orf` convention (e.g. C1orf112) and hyphenated
# read-through / antisense symbols (e.g. ANKRD13C-DT).
FORMAT_RE = re.compile(r"^[A-Z0-9][A-Za-z0-9._-]*$")


# ---------------------------------------------------------------------------
def load_hgnc(path: str | Path | None = None, refresh: bool = False) -> dict:
    """Load the HGNC complete set; download + cache if no local path given.

    Returns {"approved": set, "prev2app": {prev_symbol: approved},
             "alias": set}. Only 'Approved' records are used.
    """
    if path is not None:
        src = Path(path)
        if not src.exists():
            raise FileNotFoundError(f"HGNC file not found: {src}")
    else:
        src = HGNC_CACHE
        if refresh or not src.exists():
            print(f"[hgnc] downloading complete set -> {src}")
            urllib.request.urlretrieve(HGNC_URL, src)
        else:
            print(f"[hgnc] using cached set {src}")

    h = pd.read_csv(src, sep="\t", dtype=str, keep_default_na=False,
                    usecols=["symbol", "status", "alias_symbol", "prev_symbol"])
    h = h[h["status"] == "Approved"]

    approved = set(h["symbol"])
    prev2app: dict[str, str] = {}
    alias: set[str] = set()
    for sym, prev, ali in zip(h["symbol"], h["prev_symbol"], h["alias_symbol"]):
        for p in prev.split("|"):
            if p:
                prev2app.setdefault(p, sym)
        for a in ali.split("|"):
            if a:
                alias.add(a)
    return {"approved": approved, "prev2app": prev2app, "alias": alias}


# ---------------------------------------------------------------------------
def check_format(genes, verbose: bool = True) -> dict:
    """Gene-name formatting report (no network / no HGNC lookup)."""
    g = pd.Series(list(genes), dtype="object").astype(str)
    nonstd = sorted(g[~g.str.match(FORMAT_RE)].unique())
    report = {
        "n": len(g),
        "unique": int(g.nunique()),
        "whitespace": int((g != g.str.strip()).sum()),
        "space": int(g.str.contains(" ").sum()),
        "empty_or_na": int(((g.str.strip() == "") | (g == "NA")).sum()),
        "duplicates": int(g.duplicated().sum()),
        "nonstandard": nonstd,
        "special_chars": sorted({c for tok in g for c in tok if not c.isalnum()}),
    }
    if verbose:
        print("--- gene-name formatting check ---")
        print(f"genes: {report['n']} | unique: {report['unique']}")
        print(f"leading/trailing whitespace: {report['whitespace']}")
        print(f"embedded space: {report['space']}")
        print(f"empty or 'NA': {report['empty_or_na']}")
        print(f"duplicate symbols: {report['duplicates']}")
        print(f"non-standard-shape symbols: {len(report['nonstandard'])}"
              + (f" -> {report['nonstandard'][:15]}" if report["nonstandard"] else ""))
        print(f"special characters present: {report['special_chars']}")
    return report


# ---------------------------------------------------------------------------
def classify_hgnc(genes, hgnc: dict, verbose: bool = True) -> pd.Series:
    """Classify each gene: approved / previous_symbol / alias_symbol / not_found."""
    approved, prev2app, alias = hgnc["approved"], hgnc["prev2app"], hgnc["alias"]

    def cl(x: str) -> str:
        if x in approved:
            return "approved"
        if x in prev2app:
            return "previous_symbol"
        if x in alias:
            return "alias_symbol"
        return "not_found"

    g = pd.Series(list(genes), dtype="object").astype(str)
    st = g.map(cl)
    if verbose:
        print("--- HGNC status ---")
        vc = st.value_counts()
        for k in ["approved", "previous_symbol", "alias_symbol", "not_found"]:
            if k in vc:
                print(f"{k}: {int(vc[k])}")
        print(f"approved %: {round(100 * (st == 'approved').mean(), 1)}")
    return st


# ---------------------------------------------------------------------------
def harmonize(df: pd.DataFrame, hgnc: dict, gene_col: str = "gene",
              verbose: bool = True):
    """Rename outdated (previous) symbols to current HGNC-approved symbols.

    A previous symbol whose current name already exists among the input genes is
    a collision: the deprecated row is dropped so the existing correct row is
    kept, rather than creating a duplicate. Returns (new_df, info) where info has
    'renamed' (count) and 'dropped' (list of dropped deprecated symbols).
    """
    prev2app, approved = hgnc["prev2app"], hgnc["approved"]
    present = set(df[gene_col])

    def target(x: str) -> str:
        return prev2app[x] if (x not in approved and x in prev2app) else x

    tgt = df[gene_col].map(target)
    is_prev = tgt != df[gene_col]
    collide = is_prev & tgt.isin(present)          # target already exists → drop
    dropped = sorted(df.loc[collide, gene_col].unique())

    out = df.loc[~collide].copy()
    out[gene_col] = out[gene_col].map(target)
    n_renamed = int((is_prev & ~collide).sum())

    # Guard: two different previous symbols mapping to the same new target.
    if out[gene_col].duplicated().any():
        dups = sorted(out.loc[out[gene_col].duplicated(keep=False), gene_col].unique())
        raise ValueError(f"harmonize produced duplicate gene(s): {dups}")

    info = {"renamed": n_renamed, "dropped": dropped}
    if verbose:
        print("--- HGNC harmonization ---")
        print(f"symbols renamed: {n_renamed} | deprecated dropped: {len(dropped)}"
              + (f" ({', '.join(dropped)})" if dropped else ""))
    return out, info


# ---------------------------------------------------------------------------
def _main(argv: list[str]) -> None:
    if not argv:
        sys.exit("usage: python gene_name_utils.py <csv> [gene_col]")
    path = argv[0]
    gene_col = argv[1] if len(argv) > 1 else "gene"
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    if gene_col not in df.columns:
        sys.exit(f"column '{gene_col}' not in {path}; columns: {list(df.columns)}")
    genes = df[gene_col]
    check_format(genes)
    print()
    classify_hgnc(genes, load_hgnc())


if __name__ == "__main__":
    _main(sys.argv[1:])
