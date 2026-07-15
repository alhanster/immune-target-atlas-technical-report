#!/usr/bin/env python
"""
kNN cosine-similarity approved drug-target score — CD4+ T-cell Perturb-seq
=====================================================================
Scores every perturbed gene by how similar its transcriptome-wide knockdown
signature is to that of known approved drug targets, using a
compare-then-aggregate k-nearest-neighbour rule:

    score(gene g) = max over known approved targets t != g of
                    cosine( signature(g), signature(t) )

The reference set is the approved target list
(data/reference/approved_target_genes.txt; ~723 genes), loaded by default,
and is the ONLY set the cosine similarity is computed against. Gene symbols on
both sides are harmonized to current HGNC-approved symbols (via
Data/Gene List/gene_name_utils.py) before matching, so perturb-seq rows carrying
a deprecated *previous* symbol still match the approved anchors instead of
being silently dropped.

Rationale for MAX (not mean/centroid): the approved-target reference set is
mechanistically MULTI-MODAL (JAK/STAT, calcineurin, nucleotide synthesis,
co-stimulation, cytokine receptors ...). Max keeps the single best mechanistic
match and ignores irrelevant distances; averaging over the whole set dilutes a
real one-pathway match with unrelated positives. Empirically max beats top-k
mean and the Rocchio centroid on AP and early enrichment.

INPUT (streamed, not fully downloaded)
  s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad
  33,983 (perturbation x condition) rows x 10,282 measured genes.
  Only the contiguous, uncompressed `zscore` layer (~2.8 GB float64) is pulled
  via HTTP range requests into a disk-backed float32 memmap (~1.4 GB); the full
  16.8 GB file is never downloaded or held in RAM.

WHY zscore (= log2FC / lfcSE): variance-stabilised, so cosine similarity is not
dominated by noisy low-expression genes the way raw log_fc would be.

OUTPUT  one CSV per condition:  knn_immune_target_candidates_<COND>.csv
  columns: gene, kNN_max_score, nearest_approved_target, is_approved_target,
           n_downstream_z3, culture_condition

USAGE
  python knn_immune_target_score.py                       # Stim8hr + Stim48hr
  python knn_immune_target_score.py --conditions Rest Stim8hr Stim48hr
  python knn_immune_target_score.py --approved_targets approved_targets.txt # custom reference list
Requires: numpy, pandas, scipy, h5py  (~2 GB free disk, ~2 GB RAM)
         + Data/Gene List/gene_name_utils.py (HGNC harmonization; cached, offline)
"""
from __future__ import annotations
import argparse, io, os, sys, time, urllib.request
from pathlib import Path
import numpy as np, pandas as pd

# This script lives in "Part 4 - KNN/"; the shared gene-name utilities and the
# approved reference list live under Data/ one directory up.
HERE          = Path(__file__).resolve().parent        # src/knn
REPO          = HERE.parents[1]                         # repo root
DEFAULT_APPROVED   = REPO / "data" / "reference" / "approved_target_genes.txt"
OUTDIR        = REPO / "data" / "derived"
sys.path.insert(0, str(REPO / "src" / "shared"))
from gene_name_utils import load_hgnc, harmonize  # noqa: E402

S3 = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
      "marson2025_data/GWCD4i.DE_stats.h5ad")
MEMMAP = str(REPO / "data" / "raw" / "zscore_f32.dat")   # ~1.3 GB cache (gitignored)
NROW, NCOL = 33983, 10282


# ----------------------------------------------------------------------
# HTTP range-request file: stream arbitrary byte ranges of the remote .h5ad
# ----------------------------------------------------------------------
class HTTPRangeFile(io.RawIOBase):
    def __init__(self, url, block=4 * 1024 * 1024):
        self.url, self.block, self.pos, self.cache = url, block, 0, {}
        with urllib.request.urlopen(
                urllib.request.Request(url, method="HEAD"), timeout=60) as r:
            self.size = int(r.headers["Content-Length"])
    def _fetch(self, b):
        if b in self.cache: return self.cache[b]
        s = b * self.block; e = min(s + self.block, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={s}-{e}"})
        with urllib.request.urlopen(req, timeout=120) as r: data = r.read()
        if len(self.cache) > 256: self.cache.clear()
        self.cache[b] = data; return data
    def seek(self, off, whence=0):
        self.pos = (off if whence == 0 else
                    self.pos + off if whence == 1 else self.size + off)
        return self.pos
    def tell(self): return self.pos
    def seekable(self): return True
    def readable(self): return True
    def read(self, n=-1):
        if n < 0: n = self.size - self.pos
        end = min(self.pos + n, self.size); out = bytearray()
        for b in range(self.pos // self.block, (end - 1) // self.block + 1):
            data = self._fetch(b); s = b * self.block
            out += data[max(self.pos, s) - s: min(end, s + len(data)) - s]
        self.pos = end; return bytes(out)
    def readinto(self, b):
        d = self.read(len(b)); b[:len(d)] = d; return len(d)


def open_h5ad():
    import h5py
    return h5py.File(HTTPRangeFile(S3), "r")


def read_cat(hf, grp):
    """Read an AnnData categorical (categories+codes) or plain array."""
    import h5py
    o = hf[grp]
    if isinstance(o, h5py.Group) and "categories" in o:
        return np.asarray(o["categories"].asstr()[:])[o["codes"][:]]
    return o.asstr()[:] if o.dtype.kind in "OS" else o[:]


def stream_zscore(hf, refresh=False):
    """Stream the contiguous zscore layer into a float32 memmap (cached)."""
    lfc = hf["layers"]["zscore"]
    nrow, ncol = lfc.shape
    off = lfc.id.get_offset()
    assert lfc.chunks is None and off is not None, "expect contiguous layout"
    if refresh and os.path.exists(MEMMAP):
        os.remove(MEMMAP)
    if not (os.path.exists(MEMMAP) and os.path.getsize(MEMMAP) == nrow * ncol * 4):
        mm = np.memmap(MEMMAP, dtype=np.float32, mode="w+", shape=(nrow, ncol))
        step, t0 = 400, time.time()
        for r0 in range(0, nrow, step):
            r1 = min(r0 + step, nrow)
            b0, b1 = off + r0 * ncol * 8, off + r1 * ncol * 8 - 1
            for att in range(5):
                try:
                    req = urllib.request.Request(S3, headers={"Range": f"bytes={b0}-{b1}"})
                    with urllib.request.urlopen(req, timeout=180) as rr: buf = rr.read()
                    break
                except Exception:
                    if att == 4: raise
                    time.sleep(2 * (att + 1))
            mm[r0:r1, :] = np.frombuffer(buf, dtype="<f8").astype(np.float32).reshape(r1 - r0, ncol)
            if r0 % 8000 < step:
                print(f"[zscore] {r0}/{nrow} rows  {time.time()-t0:.0f}s", flush=True)
        mm.flush()
    return np.memmap(MEMMAP, dtype=np.float32, mode="r", shape=(nrow, ncol))


def make_canon(hgnc):
    """Return a no-drop symbol canonicaliser: previous HGNC symbol -> approved.

    Same previous->approved rule harmonize() applies internally, but element-wise
    so it can map the perturb-seq gene array WITHOUT dropping collision rows
    (which would desync `rg` from the memmap rows). Aliases are left untouched
    (ambiguous), matching harmonize()'s deliberate convention.
    """
    prev2app, approved = hgnc["prev2app"], hgnc["approved"]
    return lambda x: prev2app[x] if (x not in approved and x in prev2app) else x


# ----------------------------------------------------------------------
# Core: kNN max-score for one culture condition
# ----------------------------------------------------------------------
def score_condition(COND, mm, cond, pert_name, ont_sig, n_cells,
                    var_name, reference, canon, min_cells=50):
    name2col = {g: i for i, g in enumerate(var_name)}

    # Quality filter: this condition, significant on-target KD, enough cells.
    sel = (cond == COND) & (ont_sig.astype(bool)) & (n_cells >= min_cells)
    ridx = np.where(sel)[0]
    rg = pert_name[ridx]
    _, keep = np.unique(rg, return_index=True)          # de-duplicate gene rows
    ridx = ridx[np.sort(keep)]; rg = pert_name[ridx]

    # Signatures; zero each gene's own on-target column (drop self-knockdown).
    X = np.asarray(mm[ridx, :], dtype=np.float32)
    for r, g in enumerate(rg):
        c = name2col.get(g)
        if c is not None: X[r, c] = 0.0
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)   # cosine = dot

    S = Xn @ Xn.T
    np.fill_diagonal(S, -np.inf)                        # never self-compare
    # Harmonize perturb-seq symbols to approved HGNC for membership tests only,
    # so rows carrying a deprecated symbol still match the approved anchors.
    rg_c = np.array([canon(g) for g in rg])
    pos = np.array([g in reference for g in rg_c])      # approved-anchored positives
    pix = np.where(pos)[0]

    # compare-then-aggregate: MAX cosine to any approved target (LOO for positives)
    sc = np.empty(len(rg)); nn = np.empty(len(rg), dtype=object)
    for i in range(len(rg)):
        cols = pix[pix != i]
        j = cols[np.argmax(S[i, cols])]
        sc[i] = S[i, j]; nn[i] = rg[j]

    strength = (np.abs(X) > 3).sum(1)                   # trans-effect breadth
    return pd.DataFrame({
        "gene": rg, "kNN_max_score": sc, "nearest_approved_target": nn,
        "is_approved_target": pos, "n_downstream_z3": strength,
        "culture_condition": COND,
    }).sort_values("kNN_max_score", ascending=False).reset_index(drop=True)


def main(conditions, approved_path, refresh=False, min_cells=50):
    # Reference (anchor) set = approved targets, harmonized to current HGNC.
    # The list is already 100% approved today, so harmonize is a no-op here;
    # it future-proofs against edits to the list.
    hgnc = load_hgnc()
    canon = make_canon(hgnc)
    approved_lines = [l.strip() for l in open(approved_path or DEFAULT_APPROVED) if l.strip()]
    reference = set(harmonize(pd.DataFrame({"gene": approved_lines}), hgnc)[0]["gene"])
    print(f"[ref] {len(reference)} approved reference genes")

    hf = open_h5ad()
    cond      = read_cat(hf, "obs/culture_condition")
    pert_name = read_cat(hf, "obs/target_contrast_gene_name")
    ont_sig   = read_cat(hf, "obs/ontarget_significant")
    n_cells   = np.asarray(hf["obs/n_cells_target"][:], float)
    var_name  = np.asarray(hf["var/gene_name"].asstr()[:])
    mm = stream_zscore(hf, refresh=refresh)

    for C in conditions:
        df = score_condition(C, mm, cond, pert_name, ont_sig, n_cells,
                             var_name, reference, canon, min_cells)
        out = str(OUTDIR / f"knn_immune_target_candidates_{C}.csv")
        df.to_csv(out, index=False)
        print(f"[{C}] {len(df)} genes ({int(df.is_approved_target.sum())} approved)  ->  {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--conditions", nargs="+", default=["Stim8hr", "Stim48hr"],
                    choices=["Rest", "Stim8hr", "Stim48hr"])
    ap.add_argument("--approved_targets", default=None,
                    help="text file of reference target gene symbols (one per "
                         "line) to use as the kNN anchor / positive set; if "
                         "omitted, defaults to data/reference/"
                         "approved_target_genes.txt")
    ap.add_argument("--min-cells", type=int, default=50)
    ap.add_argument("--refresh", action="store_true", help="re-stream zscore memmap")
    a = ap.parse_args()
    main(a.conditions, a.approved_targets, refresh=a.refresh, min_cells=a.min_cells)
