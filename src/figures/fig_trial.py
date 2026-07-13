#!/usr/bin/env python3
"""
Reproduce the clinical-development validation figure (trial_fig.png).

Panel a: for each of the top 25 unlabeled PU nominations, the highest clinical
         stage of any drug against that target (Open Targets curated known
         drugs), with a red star marking targets that have a drug pursued for an
         "immune system disorder" indication.
Panel b: fraction of genes with a drug in clinical trials, top-25 nominations
         vs 50 random lower-ranked controls, for any indication and for immune
         indications, with Fisher exact p-values.

Data (two small files, regenerable from the model output + Open Targets):
  - trial_validation_top25.csv   : per-gene stage/drug/indication for the top 25
  - trial_validation_counts.json : aggregate top-vs-control counts for panel b

Usage:
  python3 make_trial_fig.py            # render from the cached data files
  python3 make_trial_fig.py --fetch    # re-pull from Open Targets, then render

The --fetch path (fetch_data below) reproduces the full validation pull: it reads
the genome-wide nomination table (nominations.tsv), selects the top 25 unlabeled
nominations and a seed-matched control set of 50 lower-ranked genes, queries the
Open Targets Platform GraphQL API (drugAndClinicalCandidates on the Target
entity) for each, and writes the two data files above. If those files are absent,
the script fetches automatically.

Styling uses the `figure-style` skill helpers when available and falls back to
inlined equivalents otherwise, so it runs standalone. Because the fallback
styling differs from the skill (font sizes, DPI), a standalone run reproduces
the same data and layout but is not pixel-identical to the manuscript figure;
run it in a kernel with the `figure-style` skill loaded for a pixel-matched
render.
"""
import os
import sys
import time
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

# ---- repo-relative paths (this script lives in src/figures/) ----
from pathlib import Path as _Path
_REPO = _Path(__file__).resolve().parents[2]
NOMINATIONS_TSV = str(_REPO / "outputs" / "nominations.tsv")   # genome-wide PU output (input to fetch)
OT_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
N_TOP = 25            # top unlabeled nominations to validate
N_CONTROL = 50        # random control genes sampled from rank > CONTROL_RANK_MIN
CONTROL_RANK_MIN = 200
SEED = 1234
IMMUNE_TA = "immune system disorder"     # Open Targets therapeutic-area name

# Open Targets returns clinical stage as e.g. 'PHASE_2', 'PHASE_1_2', 'APPROVAL'
_STAGE_FROM_OT = {"APPROVAL": ("Approved", 4), "PHASE_4": ("Phase IV", 4),
                  "PHASE_3": ("Phase III", 3), "PHASE_2_3": ("Phase II/III", 2.5),
                  "PHASE_2": ("Phase II", 2), "PHASE_1_2": ("Phase I/II", 1.5),
                  "PHASE_1": ("Phase I", 1), "EARLY_PHASE_1": ("Early Phase I", 0.5),
                  "PHASE_0": ("", 0), "PRECLINICAL": ("", 0)}

_KNOWN_DRUGS_QUERY = """
query($id: String!){
  target(ensemblId:$id){
    approvedSymbol
    drugAndClinicalCandidates{
      count
      rows{
        maxClinicalStage
        drug{ id name drugType maximumClinicalStage }
        diseases{ disease{ id name therapeuticAreas{ id name } } }
      }
    }
  }
}
"""


def _ot_known_drugs(ensembl_id, session):
    """Query Open Targets for drugs annotated against one target.

    Returns (n_drugs, max_stage_label, max_stage_num, max_drug,
             immune_stage_label, immune_stage_num, immune_drug, immune_disease).
    """
    r = session.post(OT_GRAPHQL_URL,
                     json={"query": _KNOWN_DRUGS_QUERY, "variables": {"id": ensembl_id}},
                     timeout=30)
    r.raise_for_status()
    tgt = (r.json().get("data") or {}).get("target")
    if not tgt or not tgt.get("drugAndClinicalCandidates"):
        return 0, "", -1, "", "", -1, "", ""
    dac = tgt["drugAndClinicalCandidates"]
    best = (-1, "", ""); ibest = (-1, "", "", "")      # (num, label, drug[, disease])
    for row in dac["rows"]:
        label, num = _STAGE_FROM_OT.get(row["maxClinicalStage"], ("", -1))
        drug = (row["drug"] or {}).get("name", "")
        tas, dnames = set(), []
        for d in row["diseases"]:
            dd = d.get("disease") or {}
            if dd.get("name"):
                dnames.append(dd["name"])
            for ta in (dd.get("therapeuticAreas") or []):
                tas.add(ta["name"])
        if num > best[0]:
            best = (num, label, drug)
        if IMMUNE_TA in tas and num > ibest[0]:
            ibest = (num, label, drug, dnames[0] if dnames else "")
    return (dac["count"], best[1], best[0], best[2],
            ibest[1], ibest[0], ibest[2], ibest[3])


def fetch_data():
    """Reproduce the full validation pull and write the two data files.

    Requires `requests` and the genome-wide nomination table NOMINATIONS_TSV.
    """
    import requests
    nom = pd.read_csv(NOMINATIONS_TSV, sep="\t")

    # top N unlabeled nominations (not approved immune-drug targets)
    top = (nom[nom.label_immune == 0].sort_values("rank").head(N_TOP)
           [["rank", "gene", "gene_id", "pu_score"]].reset_index(drop=True))
    # matched controls: random unlabeled genes ranked outside the top CONTROL_RANK_MIN
    pool = nom[(nom.label_immune == 0) & (nom["rank"] > CONTROL_RANK_MIN)]
    ctrl = pool.sample(N_CONTROL, random_state=SEED)[["rank", "gene", "gene_id"]]

    session = requests.Session()
    rows = []
    for _, r in top.iterrows():
        (n, ms, msn, md, ims, imsn, imd, imdis) = _ot_known_drugs(r["gene_id"], session)
        rows.append(dict(pu_rank=int(r["rank"]), gene=r["gene"], n_known_drugs=n,
                         max_clinical_stage=ms, example_drug=md,
                         immune_max_stage=ims, immune_example_drug=imd,
                         immune_indication=imdis))
        time.sleep(0.15)                       # be polite to the shared API
    top_df = pd.DataFrame(rows)
    top_df.to_csv(TOP_CSV, index=False)

    # aggregate control counts (only the tallies are needed for panel b)
    c_any = c_imm = 0
    for _, r in ctrl.iterrows():
        (n, ms, msn, md, ims, imsn, imd, imdis) = _ot_known_drugs(r["gene_id"], session)
        c_any += int(n > 0)
        c_imm += int(imsn >= 0)
        time.sleep(0.15)

    counts = {"top_n": int(len(top_df)),
              "top_any_drug": int((top_df.n_known_drugs > 0).sum()),
              "top_immune": int(top_df.immune_max_stage.astype(str).str.strip().ne("").sum()),
              "ctrl_n": int(len(ctrl)),
              "ctrl_any_drug": int(c_any),
              "ctrl_immune": int(c_imm)}
    json.dump(counts, open(COUNTS_JSON, "w"), indent=2)
    print(f"fetched Open Targets data -> {TOP_CSV}, {COUNTS_JSON}")

# --- figure styling (skill helpers, with standalone fallback) ---
try:
    from figure_style import apply_figure_style  # noqa
except ImportError:
    def apply_figure_style(sizes=(9, 8, 7)):
        base, mid, small = sizes
        plt.rcParams.update({
            "font.size": base, "axes.titlesize": base, "axes.labelsize": mid,
            "xtick.labelsize": small, "ytick.labelsize": small,
            "axes.spines.top": False, "axes.spines.right": False,
            "figure.dpi": 110,
        })

TOP_CSV = str(_REPO / "data" / "derived" / "trial_validation_top25.csv")
COUNTS_JSON = str(_REPO / "data" / "derived" / "trial_validation_counts.json")
OUT_PNG = str(_REPO / "plots" / "trial_fig.png")

# ordinal phase -> numeric position on the panel-a axis
STAGE_NUM = {"Approved": 4, "Phase IV": 4, "Phase III": 3, "Phase II/III": 2.5,
             "Phase II": 2, "Phase I/II": 1.5, "Phase I": 1, "Early Phase I": 0.5, "": -1}
# single-hue ordinal ramp: deeper blue = later stage
STAGE_COL = {"Approved": "#08519c", "Phase IV": "#08519c", "Phase III": "#2a7ab0",
             "Phase II/III": "#4a97c9", "Phase II": "#74add1", "Phase I/II": "#a6bddb",
             "Phase I": "#c6dbef", "": "#e3e3e3"}


def main():
    # fetch from Open Targets if asked, or if the cached data files are missing
    need_fetch = "--fetch" in sys.argv or not (os.path.exists(TOP_CSV) and os.path.exists(COUNTS_JSON))
    if need_fetch:
        fetch_data()

    apply_figure_style(sizes=(9, 8, 7))

    top = pd.read_csv(TOP_CSV).fillna("")
    top = top.sort_values("pu_rank").reset_index(drop=True)
    top["max_num"] = top["max_clinical_stage"].map(lambda s: STAGE_NUM.get(str(s), -1))
    top["is_immune"] = top["immune_max_stage"].astype(str).str.strip().ne("") \
        & top["immune_max_stage"].astype(str).str.strip().ne("nan")

    counts = json.load(open(COUNTS_JSON))
    tn, cn = counts["top_n"], counts["ctrl_n"]
    t_any, c_any = counts["top_any_drug"], counts["ctrl_any_drug"]
    t_imm, c_imm = counts["top_immune"], counts["ctrl_immune"]
    _, p_any = fisher_exact([[t_any, tn - t_any], [c_any, cn - c_any]])
    _, p_imm = fisher_exact([[t_imm, tn - t_imm], [c_imm, cn - c_imm]])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.8, 5.2),
                                   gridspec_kw={"width_ratios": [2.1, 1]})

    # ---------- Panel a: per-gene highest clinical stage ----------
    y = np.arange(len(top))[::-1]
    for yi, (_, r) in zip(y, top.iterrows()):
        if r["max_num"] < 0:
            axA.barh(yi, 0.10, color="#f0f0f0", edgecolor="#cfcfcf",
                     height=0.68, linewidth=0.5)          # stub = no known drug
        else:
            axA.barh(yi, r["max_num"], color=STAGE_COL.get(r["max_clinical_stage"], "#e3e3e3"),
                     height=0.68, edgecolor="white", linewidth=0.4)
        if r["is_immune"]:
            axA.scatter(max(r["max_num"], 0) + 0.16, yi, marker="*", s=64,
                        color="#d1495b", zorder=5)
    axA.set_yticks(y)
    axA.set_yticklabels(top["gene"], fontstyle="italic", fontsize=6.6)
    axA.set_xticks([0, 1, 2, 3, 4])
    axA.set_xticklabels(["none", "Ph I", "Ph II", "Ph III", "Appr"])
    axA.set_xlabel("Highest clinical stage of a drug against the target", labelpad=8)
    axA.set_xlim(0, 4.7)
    axA.set_ylim(-0.8, len(top) - 0.2)
    axA.set_title("Top 25 novel nominations: existing drug programs", loc="left", fontsize=9)
    # in-panel legend, placed in empty lower-mid whitespace
    axA.scatter([2.55], [4.0], marker="*", s=64, color="#d1495b", zorder=5)
    axA.text(2.75, 4.0, "= in trials for an\n   immune indication",
             fontsize=6.6, va="center", ha="left")

    # ---------- Panel b: enrichment vs control ----------
    cats = ["Any drug\nin trials", "Immune-indication\ntrial drug"]
    top_pct = [100 * t_any / tn, 100 * t_imm / tn]
    ctrl_pct = [100 * c_any / cn, 100 * c_imm / cn]
    x = np.arange(len(cats))
    w = 0.36
    axB.bar(x - w / 2, top_pct, w, color="#2a7ab0", label=f"Top {tn} nominations")
    axB.bar(x + w / 2, ctrl_pct, w, color="#b0b0b0", label=f"Random controls (n={cn})")
    for xi, tp, cp in zip(x, top_pct, ctrl_pct):
        axB.text(xi - w / 2, tp + 0.8, f"{tp:.0f}%", ha="center",
                 fontsize=7.5, fontweight="bold")
        axB.text(xi + w / 2, cp + 0.8, f"{cp:.0f}%", ha="center",
                 fontsize=7.5, color="#666")
    axB.set_xticks(x)
    axB.set_xticklabels(cats, fontsize=7)
    axB.set_ylabel("% of genes")
    axB.set_ylim(0, 30)
    axB.set_title("Nominations enriched for\nactive drug programs", loc="left", fontsize=9)
    axB.legend(frameon=False, fontsize=6.2, loc="upper right")
    axB.annotate(f"p={p_any:.3f}", xy=(0, 26.5), fontsize=6.4, ha="center", color="#08519c")
    axB.annotate(f"p={p_imm:.3f}", xy=(1, 24.5), fontsize=6.4, ha="center", color="#08519c")

    fig.tight_layout(pad=1.0)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT_PNG}  (any-drug p={p_any:.3f}, immune p={p_imm:.3f})")


if __name__ == "__main__":
    main()
