#!/usr/bin/env python
"""
Genome-wide regulator-burden ("core gene") scoring pipeline
===========================================================
Scores every gene as a candidate *core gene* for lymphocyte count by
integrating a CD4+ T-cell genome-scale Perturb-seq DE matrix with UK Biobank
loss-of-function burden estimates. No gene panel is involved -- scoring is
genome-wide; any downstream gene-set selection (e.g. IEI) is done by the
consumer (see make_regulator_burden_figure.R, Panel b).

Method (per candidate core gene j, per culture condition), a Python port of
the authors' reference R script
  GWT_perturbseq_analysis_2025/src/8_lymphocyte_counts_LoF/
  Regulator_burden_correlation_GWT.R :

    scale(gamma_x) ~ a * scale(beta_{x->j}) + b * shet_x + intercept   (x != j)

  beta_{x->j} = log2FC of gene j after CRISPRi knockdown of regulator x
                (from the Perturb-seq DE matrix)
  gamma_x     = GeneBayes-denoised LoF burden effect of x on lymphocyte
                count (UKB/Backman). This is the burden EFFECT SIZE on the
                trait -- NOT the raw lymphocyte count, and NOT perturb-seq.
  shet_x      = selective-constraint nuisance covariate
  Self-effect (x == j) is masked.

The per-gene score for j is coefficient `a` (`coef_beta`): sign = direction
(positive => up-regulation of j promotes the trait), p-value = significance.

DELIVERABLES (written to the working directory)
  regulator_burden_scores_{cond}.csv     genome-wide scores per condition
  regulator_burden_scores_all.csv        genome-wide scores, all conditions
                                         stacked (main deliverable)

The manuscript figure (signed-QQ + core-gene bar chart) is drawn separately by
make_regulator_burden_figure.R, which reads these CSVs -- this script does
not produce figures.

DATA SOURCES
  1. Analysis repo  https://github.com/emdann/GWT_perturbseq_analysis_2025
     ships: gamma (LoF burden), S_het, the Ensembl<->symbol map, and the
     authors' own core-gene lists (validation).
  2. CZI Virtual Cell Models bucket (the DE matrix beta = log_fc layer):
       s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad
     Public anonymous read, 16.8 GB. Only the contiguous log_fc byte range
     (~2.8 GB) is streamed via HTTP range requests and cached to disk as a
     float32 memmap -- the full file is never downloaded or held in RAM.

RUNTIME  (dominated by the one-time beta stream)
  --self-test   seconds,  no network/data   -> validates the scoring math
  (default)     ~3 min first run / ~40 s cached  (genome-wide scoring)
  The beta memmap (logfc_f32.dat, ~1.4 GB) is cached; re-runs skip the
  download. Delete it or pass --refresh to re-stream.

REQUIREMENTS   python -m pip install -r requirements.txt   (+ git on PATH)
  ~2 GB free disk for the cached matrix; ~2 GB RAM.
  Network: github.com and genome-scale-tcell-perturb-seq.s3.amazonaws.com

USAGE
  python regulator_burden_pipeline.py --self-test
  python regulator_burden_pipeline.py            # full, genome-wide
"""
from __future__ import annotations
import argparse, io, os, time, urllib.request
import numpy as np
import pandas as pd
from scipy import stats as st

# Small metadata inputs are vendored into inputs/ (see inputs/SOURCE.md);
# they originate from the MIT-licensed repo emdann/GWT_perturbseq_analysis_2025
# and its upstream data sources. Only the large DE matrix (beta) is fetched at
# runtime from the CZI bucket.
S3_H5AD  = ("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
            "marson2025_data/GWCD4i.DE_stats.h5ad")
# Paths anchored to this script's location so it runs from any cwd. The vendored
# vendored metadata inputs (and the cached beta memmap) live in this stage's inputs/.
HERE     = os.path.dirname(os.path.abspath(__file__))   # src/regulator_burden
REPO     = os.path.dirname(os.path.dirname(HERE))        # repo root
INPUT    = os.path.join(HERE, "inputs")
RESULTS  = os.path.join(INPUT, "core_genes_reference")
MEMMAP   = os.path.join(INPUT, "logfc_f32.dat")          # ~1.3 GB cache (gitignored)
OUTDIR   = os.path.join(REPO, "data", "derived")         # deliverables -> data/derived/
CONDS    = ["Rest", "Stim8hr", "Stim48hr"]


# ----------------------------------------------------------------------
# 1. Verify vendored metadata inputs are present
# ----------------------------------------------------------------------
def check_inputs():
    needed = [
        f"{INPUT}/Backman_LymphocyteCount_fullFeatures.per_gene_estimates.tsv",
        f"{INPUT}/shet_10bins.txt",
        f"{INPUT}/gencode_v41_gname_gid_ALL_sorted",
    ]
    missing = [f for f in needed if not os.path.exists(f)]
    if missing:
        raise SystemExit(
            "Missing vendored input files:\n  " + "\n  ".join(missing) +
            f"\nExpected under {INPUT}/ (see inputs/SOURCE.md for provenance).")
    print(f"[inputs] all vendored metadata present under {INPUT}/")


# ----------------------------------------------------------------------
# 2. Stream the log_fc layer from S3 into a disk-backed float32 memmap
#    Reads the .h5ad header via HTTP range requests to locate the
#    contiguous dataset, then streams only that byte range.
# ----------------------------------------------------------------------
class HTTPRangeFile(io.RawIOBase):
    """Seekable read-only file over HTTP range requests (block-cached)."""
    def __init__(self, url, block=4 * 1024 * 1024):
        self.url, self.block, self.pos, self.cache = url, block, 0, {}
        with urllib.request.urlopen(
                urllib.request.Request(url, method="HEAD"), timeout=60) as r:
            self.size = int(r.headers["Content-Length"])
    def _fetch(self, b):
        if b in self.cache: return self.cache[b]
        s = b * self.block; e = min(s + self.block, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={s}-{e}"})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        if len(self.cache) > 256: self.cache.clear()
        self.cache[b] = data
        return data
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
        self.pos = end
        return bytes(out)
    def readinto(self, b):
        d = self.read(len(b)); b[:len(d)] = d; return len(d)


def load_de_matrix():
    """Return (memmap[obs,var] float32 log_fc, obs_condition, obs_pert_name,
    var_gene_ids, var_gene_name)."""
    import h5py
    hf = h5py.File(HTTPRangeFile(S3_H5AD), "r")
    lfc = hf["layers"]["log_fc"]
    nrow, ncol = lfc.shape
    assert lfc.dtype.byteorder in "<=|", "expect little-endian float64"

    def read_cat(grp):                      # AnnData categorical or array
        o = hf[grp]
        if isinstance(o, h5py.Group) and "categories" in o:
            return o["categories"].asstr()[:][o["codes"][:]]
        return o.asstr()[:] if o.dtype.kind in "OS" else o[:]

    cond      = read_cat("obs/culture_condition")
    pert_name = read_cat("obs/target_contrast_gene_name")
    var_ids   = hf["var"]["gene_ids"].asstr()[:]
    var_name  = hf["var"]["gene_name"].asstr()[:]

    if not (os.path.exists(MEMMAP) and
            os.path.getsize(MEMMAP) == nrow * ncol * 4):
        off = lfc.id.get_offset()           # contiguous, uncompressed dataset
        assert off is not None and lfc.chunks is None, "expect contiguous layout"
        assert lfc.id.get_storage_size() == nrow * ncol * 8
        mm = np.memmap(MEMMAP, dtype=np.float32, mode="w+", shape=(nrow, ncol))
        step, t0 = 400, time.time()
        for r0 in range(0, nrow, step):
            r1 = min(r0 + step, nrow)
            b0, b1 = off + r0 * ncol * 8, off + r1 * ncol * 8 - 1
            req = urllib.request.Request(S3_H5AD, headers={"Range": f"bytes={b0}-{b1}"})
            for att in range(5):
                try:
                    with urllib.request.urlopen(req, timeout=180) as rr:
                        buf = rr.read()
                    break
                except Exception:
                    if att == 4: raise
                    time.sleep(2 * (att + 1))
            mm[r0:r1, :] = np.frombuffer(buf, dtype="<f8").astype(np.float32).reshape(r1 - r0, ncol)
            if r0 % 8000 < step:
                print(f"[beta] {r0}/{nrow} rows  {time.time()-t0:.0f}s", flush=True)
        mm.flush()
    mm = np.memmap(MEMMAP, dtype=np.float32, mode="r", shape=(nrow, ncol))
    print(f"[beta] log_fc matrix ready {mm.shape}  range [{mm.min():.1f},{mm.max():.1f}]")
    return mm, cond, pert_name, var_ids, var_name


# ----------------------------------------------------------------------
# 3. Regulator-burden scoring
# ----------------------------------------------------------------------
def load_gamma_shet():
    gamma = pd.read_csv(f"{INPUT}/Backman_LymphocyteCount_fullFeatures.per_gene_estimates.tsv", sep="\t")
    shet  = pd.read_csv(f"{INPUT}/shet_10bins.txt", sep="\t")
    gmap  = pd.read_csv(f"{INPUT}/gencode_v41_gname_gid_ALL_sorted", sep="\t",
                        header=None, names=["ensg", "symbol"])
    g = gamma.set_index("ensg")["post_mean"].astype(float)
    g = g[~g.index.duplicated()]
    fin = g[np.isfinite(g)]; g = g.clip(fin.min(), fin.max())   # clip +/-Inf
    s = shet.set_index("ensg")["shet"].astype(float)
    s = s[~s.index.duplicated()]
    return g, s, gmap


def score_condition(COND, mm, cond, pert_ensg, var_ids, var_name,
                    gam_row, shet_row, col_block=500, min_reg=100):
    ridx = np.where(cond == COND)[0]
    gam, sh, pe = gam_row[ridx], shet_row[ridx], pert_ensg[ridx]
    base_ok = np.isfinite(gam) & np.isfinite(sh)
    ncol = mm.shape[1]; rows = []
    for c0 in range(0, ncol, col_block):
        c1 = min(c0 + col_block, ncol)
        block = np.asarray(mm[ridx.min():ridx.max() + 1, c0:c1])[ridx - ridx.min(), :]
        for jc in range(c1 - c0):
            j = c0 + jc
            beta = block[:, jc]
            ok = base_ok & np.isfinite(beta) & (pe != var_ids[j])   # mask self
            n = int(ok.sum())
            if n < min_reg: continue
            b, yv, sv = beta[ok], gam[ok], sh[ok]
            bz = (b - b.mean()) / b.std(ddof=1) if b.std(ddof=1) > 0 else b * 0
            yz = (yv - yv.mean()) / yv.std(ddof=1) if yv.std(ddof=1) > 0 else yv * 0
            X = np.column_stack([bz, sv, np.ones(n)])
            coef, *_ = np.linalg.lstsq(X, yz, rcond=None)
            resid = yz - X @ coef; dof = n - 3
            se = np.sqrt((resid @ resid) / dof * np.linalg.inv(X.T @ X)[0, 0])
            p = 2 * st.t.sf(abs(coef[0] / se), dof)
            rows.append((var_ids[j], var_name[j], n, coef[0], se, p))
    r = pd.DataFrame(rows, columns=["ensg", "gene", "n_reg", "coef_beta", "se_beta", "p_beta"])
    r["signed_log10p"] = np.sign(r.coef_beta) * -np.log10(r.p_beta.clip(lower=1e-300))
    r["fdr"] = st.false_discovery_control(r.p_beta.clip(1e-300, 1))
    return r.sort_values("signed_log10p", ascending=False).reset_index(drop=True)


# ----------------------------------------------------------------------
# 4. Validation against the authors' published core-gene lists
# ----------------------------------------------------------------------
def validate(scores):
    for cd in ["Stim8hr", "Stim48hr"]:
        for dr, sub in [("positive", scores[cd].head(50)),
                        ("negative", scores[cd].tail(50))]:
            f = f"{RESULTS}/core_genes_{cd}_{dr}.txt"
            if not os.path.exists(f): continue
            auth = set(pd.read_csv(f, header=None)[0])
            print(f"[validate] {cd} {dr}: top-50 overlap "
                  f"{len(auth & set(sub.gene))}/50")


# ----------------------------------------------------------------------
# Self-test: validate the scoring math on synthetic data (no network/data).
# Plants a signal gene (beta correlated with gamma) and a null gene; the
# regression must recover the signal (p tiny, positive coef) and reject the
# null (p > 0.05).
# ----------------------------------------------------------------------
def self_test():
    rng = np.random.default_rng(0)
    n = 3000
    pe   = np.array([f"ENSG{i:08d}" for i in range(n)])
    gam  = rng.normal(0, 0.1, n)
    shet = rng.uniform(0, 0.1, n)
    beta_sig  = 0.8 * (gam - gam.mean()) + rng.normal(0, 0.05, n)   # coupled
    beta_null = rng.normal(0, 0.3, n)                              # independent
    mm   = np.column_stack([beta_sig, beta_null]).astype(np.float32)
    cond = np.array(["Rest"] * n)
    var_ids  = np.array(["GENE_SIGNAL", "GENE_NULL"])
    var_name = var_ids.copy()
    r = score_condition("Rest", mm, cond, pe, var_ids, var_name, gam, shet,
                        min_reg=50).set_index("gene")
    print(r[["n_reg", "coef_beta", "p_beta"]].round(4).to_string())
    assert r.loc["GENE_SIGNAL", "p_beta"] < 1e-10 and r.loc["GENE_SIGNAL", "coef_beta"] > 0
    assert r.loc["GENE_NULL", "p_beta"] > 0.05
    print("self-test PASSED")


# ----------------------------------------------------------------------
def main(refresh=False):
    check_inputs()

    if refresh and os.path.exists(MEMMAP):
        os.remove(MEMMAP)
    mm, cond, pert_name, var_ids, var_name = load_de_matrix()
    gamma, shet, gmap = load_gamma_shet()
    name2ensg = dict(zip(gmap.symbol, gmap.ensg))
    pert_ensg = np.array([name2ensg.get(n, n) for n in pert_name])
    gam_row, shet_row = gamma.reindex(pert_ensg).values, shet.reindex(pert_ensg).values

    scores = {}
    for C in CONDS:
        scores[C] = score_condition(C, mm, cond, pert_ensg, var_ids, var_name,
                                     gam_row, shet_row)
        scores[C].to_csv(os.path.join(OUTDIR, f"regulator_burden_scores_{C}.csv"), index=False)
        print(f"[score] {C}: {len(scores[C])} genes, FDR<0.05 = {(scores[C].fdr<0.05).sum()}")
    validate(scores)                         # authors' lists are genome-wide top-50

    frames = []
    for C in CONDS:
        s = scores[C].copy(); s["condition"] = C
        s["rank"] = s.signed_log10p.rank(ascending=False).astype(int)
        frames.append(s)
    all_sc = pd.concat(frames, ignore_index=True)   # genome-wide, all conditions
    all_sc.to_csv(os.path.join(OUTDIR, "regulator_burden_scores_all.csv"), index=False)
    print(f"[all] {all_sc.gene.nunique()} genes scored; "
          f"{(all_sc.p_beta<0.05).sum()} gene-condition tests p<0.05")

    print("\nDONE. Deliverables written to the working directory.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--self-test", action="store_true",
                    help="validate scoring math on synthetic data (no download)")
    ap.add_argument("--refresh", action="store_true",
                    help="re-stream the beta matrix even if cached")
    a = ap.parse_args()
    if a.self_test:
        self_test()
    else:
        main(refresh=a.refresh)
