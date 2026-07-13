#!/usr/bin/env python
"""
Signature PCA export for the kNN FDA-target figure — CD4+ T-cell Perturb-seq
============================================================================
Figure-support helper for make_knn_candidates_figure.R. Rebuilds the per-gene
normalized knockdown-signature matrix `Xn` (exactly as knn_immune_target_score.py
does), reduces it to the top-50 principal components, and writes a compact CSV the
R script embeds with UMAP.

WHY a separate PCs file (not raw Xn): Xn is ~7,100 x 10,282 float32 (~290 MB per
condition); its top-50 PCs preserve the global structure UMAP needs in ~3 MB, and
keep the scoring script's "two CSVs only" output contract untouched.

INPUT   (all reused, nothing re-downloaded)
  - the cached zscore memmap (zscore_f32.dat) via knn_immune_target_score.py
  - knn_immune_target_candidates_<COND>.csv  (for is_fda_target / nearest_fda_target)

OUTPUT  one CSV per condition:  knn_signature_pcs_<COND>.csv
  columns: gene, is_fda_target, nearest_fda_target, PC1 ... PC50

USAGE
  python embed_signatures.py                       # Stim8hr + Stim48hr
  python embed_signatures.py --conditions Rest Stim8hr Stim48hr
Requires: numpy, pandas, scikit-learn, h5py  + the scorer's streaming helpers.
"""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd
from sklearn.decomposition import PCA

# Reuse the scorer's streaming helpers so the two stay in lockstep and nothing is
# re-downloaded (the memmap is cached). FDA membership / nearest-target labels are
# merged from the candidate CSV, so no HGNC/reference machinery is needed here.
from knn_immune_target_score import open_h5ad, read_cat, stream_zscore

from pathlib import Path
DERIVED = Path(__file__).resolve().parents[2] / "data" / "derived"

N_PCS = 50


def build_Xn(COND, mm, cond, pert_name, ont_sig, n_cells, var_name, min_cells=50):
    """Rebuild the L2-normalized signature matrix `Xn` and its gene order.

    MUST mirror score_condition() in knn_immune_target_score.py exactly (same QC
    filter, de-duplication, self-column zeroing, normalization) so the rows here
    correspond 1:1 to the genes in knn_immune_target_candidates_<COND>.csv.
    """
    name2col = {g: i for i, g in enumerate(var_name)}
    sel = (cond == COND) & (ont_sig.astype(bool)) & (n_cells >= min_cells)
    ridx = np.where(sel)[0]
    rg = pert_name[ridx]
    _, keep = np.unique(rg, return_index=True)          # de-duplicate gene rows
    ridx = ridx[np.sort(keep)]; rg = pert_name[ridx]

    X = np.asarray(mm[ridx, :], dtype=np.float32)
    for r, g in enumerate(rg):                          # zero each gene's own KD col
        c = name2col.get(g)
        if c is not None: X[r, c] = 0.0
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    return rg, Xn


def main(conditions, min_cells=50):
    hf = open_h5ad()
    cond      = read_cat(hf, "obs/culture_condition")
    pert_name = read_cat(hf, "obs/target_contrast_gene_name")
    ont_sig   = read_cat(hf, "obs/ontarget_significant")
    n_cells   = np.asarray(hf["obs/n_cells_target"][:], float)
    var_name  = np.asarray(hf["var/gene_name"].asstr()[:])
    mm = stream_zscore(hf)

    for C in conditions:
        rg, Xn = build_Xn(C, mm, cond, pert_name, ont_sig, n_cells, var_name, min_cells)
        pca = PCA(n_components=N_PCS, random_state=0)
        pcs = pca.fit_transform(Xn)
        df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(N_PCS)])
        df.insert(0, "gene", rg)

        # Attach the figure's color keys from the already-written candidate CSV.
        cand = pd.read_csv(DERIVED / f"knn_immune_target_candidates_{C}.csv",
                           usecols=["gene", "is_fda_target", "nearest_fda_target"])
        df = df.merge(cand, on="gene", how="left")
        df = df[["gene", "is_fda_target", "nearest_fda_target"]
                + [f"PC{i+1}" for i in range(N_PCS)]]

        out = str(DERIVED / f"knn_signature_pcs_{C}.csv")
        df.to_csv(out, index=False)
        ev = pca.explained_variance_ratio_.sum()
        print(f"[{C}] {len(df)} genes x {N_PCS} PCs ({ev:.1%} var)  ->  {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--conditions", nargs="+", default=["Stim8hr", "Stim48hr"],
                    choices=["Rest", "Stim8hr", "Stim48hr"])
    ap.add_argument("--min-cells", type=int, default=50)
    a = ap.parse_args()
    main(a.conditions, min_cells=a.min_cells)
