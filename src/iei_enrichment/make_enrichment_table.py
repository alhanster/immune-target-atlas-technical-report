"""
Build the IEI x approved-immune-drug-target enrichment table across background
universes, from the four source gene lists.

Reproduces the summary table (the uploaded photo):
    Background universe | N | IEI target rate | non-IEI rate | Odds ratio (95% CI) | Fisher p
and the readable companion CSV with counted fractions.

Source files (gene symbols):
    IEI_gene_list.csv                 505 unique inborn-errors-of-immunity genes
    approved_target_genes.txt         723 immune-drug target genes (OT-approved)
    druggable_genome_gene_list.xlsx   ~4,479 druggable-genome genes (Finan et al. 2017)
    rna_immune_cell.tsv               HPA RNA-seq, 20,162 genes x 19 sorted immune
                                      cell types (long format; nTPM column)

The "immune-expressed" universe is DERIVED, not downloaded: a gene is immune-
expressed at threshold t if its max nTPM across the 19 cell types is >= t.
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# ---- repo-relative file paths (this script lives in src/iei_enrichment/) ----
from pathlib import Path
_REPO     = Path(__file__).resolve().parents[2]
IEI_CSV   = _REPO / "data" / "reference" / "IEI_gene_list.csv"
APPROVED_TXT   = _REPO / "data" / "reference" / "approved_target_genes.txt"
DRUG_XLSX = _REPO / "data" / "reference" / "druggable_genome_gene_list.xlsx"
# HPA immune-cell RNA reduced to per-gene max nTPM (committed subset; see DATA.md).
HPA_TSV   = _REPO / "data" / "reference" / "hpa_immune_max_ntpm.tsv"
_DERIVED  = _REPO / "data" / "derived"


def load_lists():
    iei = set(pd.read_csv(IEI_CSV)['Gene'].dropna().astype(str).str.strip())
    with open(APPROVED_TXT) as f:
        approved_targets = set(l.strip() for l in f if l.strip())
    drug = set(pd.read_excel(DRUG_XLSX)['Gene'].dropna().astype(str).str.strip())
    return iei, approved_targets, drug


def immune_expressed_sets(thresholds=(1, 5, 10)):
    """Return {threshold: set(symbols)} and the full HPA symbol universe."""
    df = pd.read_csv(HPA_TSV, sep="\t")                      # Gene, Gene name, Immune cell, TPM, pTPM, nTPM
    g = df.groupby(['Gene', 'Gene name'])['nTPM'].max().reset_index()
    g.columns = ['ensembl', 'symbol', 'max_nTPM']
    sets = {t: set(g.loc[g['max_nTPM'] >= t, 'symbol'].astype(str)) for t in thresholds}
    return sets, set(g['symbol'].dropna().astype(str)), g


def enrich(iei, target, universe, label):
    """2x2 enrichment of target-status among IEI genes, restricted to `universe`."""
    U = universe
    iei_u = iei & U
    tgt_u = target & U
    a = len(iei_u & tgt_u)          # IEI & target
    b = len(iei_u) - a              # IEI & not target
    c = len(tgt_u) - a              # non-IEI & target
    d = len(U) - a - b - c          # neither
    OR, p = fisher_exact([[a, b], [c, d]], alternative='greater')
    rr = (a / len(iei_u)) / (c / (len(U) - len(iei_u))) if c > 0 else np.inf
    if min(a, b, c, d) > 0:         # Woolf 95% CI on the odds ratio
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        lo, hi = np.exp(np.log(OR) - 1.96*se), np.exp(np.log(OR) + 1.96*se)
    else:
        lo, hi = np.nan, np.nan
    return dict(universe=label, N=len(U), a=a, b=b, c=c, d=d,
                p_tgt_iei=a/len(iei_u), p_tgt_non=c/(len(U)-len(iei_u)),
                RR=rr, OR=OR, OR_lo=lo, OR_hi=hi, p=p)


def build_enrichment_table():
    iei, approved_targets, drug = load_lists()
    imm, hpa_all, _ = immune_expressed_sets()

    all_pc = hpa_all | iei | approved_targets | drug        # whole-genome proxy (~20,283 symbols)
    universes = [
        ("All genes (HPA-20k proxy)",     all_pc),
        ("Immune-expressed nTPM>=1",      imm[1]),
        ("Immune-expressed nTPM>=5",      imm[5]),
        ("Druggable genome",              drug),
        ("Druggable ∩ immune-expr(>=1)",  drug & imm[1]),
    ]
    res = pd.DataFrame([enrich(iei, approved_targets, U, lab) for lab, U in universes])

    # readable companion: every rate as a counted fraction; fold from FULL-precision
    # ratio rounded ONCE to 1 dp (avoid double-rounding, e.g. 2.35 -> 2.4)
    def fold(r):
        return f"{(r.p_tgt_iei / r.p_tgt_non):.1f}× (OR {r.OR:.1f}, {r.OR_lo:.1f}–{r.OR_hi:.1f})"
    disp = ['All genes (~20,000)', 'Immune-expressed (HPA, nTPM≥1)',
            'Immune-expressed (HPA, nTPM≥5)', 'Druggable genome (Finan et al. 2017)',
            'Druggable ∩ immune-expressed']
    readable = pd.DataFrame({
        'Comparison group': disp,
        'IEI genes that are targets':   [f"{r.a} / {r.a+r.b} = {r.p_tgt_iei:.1%}" for r in res.itertuples()],
        'Other genes that are targets': [f"{r.c} / {r.c+r.d} = {r.p_tgt_non:.1%}" for r in res.itertuples()],
        'Fold higher (odds ratio, 95% CI)': [fold(r) for r in res.itertuples()],
        'Fisher p': [f"{r.p:.0e}" for r in res.itertuples()],
    })
    return res, readable


if __name__ == "__main__":
    res, readable = build_enrichment_table()
    res.to_csv(_DERIVED / "iei_enrichment_by_universe.csv", index=False)     # full numeric
    readable.to_csv(_DERIVED / "iei_enrichment_readable.csv", index=False)   # readable
    print(readable.to_string(index=False))
