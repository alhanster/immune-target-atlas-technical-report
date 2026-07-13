"""
Two-panel manuscript figure: enrichment of IEI-associated genes among FDA
immune-drug targets, across five background universes.

  Panel A  rate bars       — % of IEI vs non-IEI genes that are targets,
                             in three headline universes, with fold-brackets
  Panel B  odds-ratio forest — OR + 95% CI across all five universes

Inputs
------
Expects a DataFrame `res` with one row per background universe and columns:
    universe, N, a, b, c, d, p_tgt_iei, p_tgt_non, RR, OR, OR_lo, OR_hi, p
where for each universe (restricted to that gene set):
    a = # IEI genes that are FDA immune-drug targets
    b = # IEI genes that are not targets
    c = # non-IEI genes that are targets
    d = # non-IEI genes that are not targets
    p_tgt_iei = a/(a+b);  p_tgt_non = c/(c+d)
    OR, OR_lo, OR_hi = odds ratio and Woolf 95% CI;  p = Fisher exact (two-sided)

`res` is produced by build_enrichment_table() in make_enrichment_table.py.

Style: uses the `figure-style` skill helpers (apply_figure_style, set_frame,
panel_letter). If that skill is not loaded, replace those calls with the
plain-matplotlib fallbacks noted inline.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# --- figure-style skill helpers, with standalone fallbacks -------------------
# In-session these come from the `figure-style` skill kernel plugin. When run as
# a plain script (no skill loaded) the fallbacks below reproduce the same look.
try:
    apply_figure_style  # noqa: F821  (injected by the skill)
except NameError:
    def apply_figure_style(sizes=(9, 8, 7)):
        base, mid, small = sizes
        plt.rcParams.update({
            'font.size': base, 'axes.titlesize': base, 'axes.labelsize': base,
            'legend.fontsize': mid, 'xtick.labelsize': small, 'ytick.labelsize': small,
            'axes.spines.top': False, 'axes.spines.right': False,
            'font.family': 'sans-serif', 'figure.dpi': 100,
        })

    def set_frame(ax):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    def panel_letter(ax, letter):
        ax.text(-0.10, 1.06, letter, transform=ax.transAxes,
                fontsize=13, fontweight='bold', va='top', ha='right')

apply_figure_style(sizes=(9, 8, 7))


def make_figure(res, outfile="iei_enrichment_manuscript.png"):
    sub = res.set_index("universe")

    order5 = ['All genes (HPA-20k proxy)', 'Immune-expressed nTPM>=1',
              'Immune-expressed nTPM>=5', 'Druggable genome',
              'Druggable ∩ immune-expr(>=1)']
    disp5  = ['All genes (~20,000)', 'Immune-expressed (HPA, nTPM≥1)',
              'Immune-expressed (HPA, nTPM≥5)', 'Druggable genome (Finan 2017)',
              'Druggable ∩ immune-expressed']
    order3 = ['All genes (HPA-20k proxy)', 'Immune-expressed nTPM>=1', 'Druggable genome']
    disp3  = ['All genes\n(~20,000)', 'Immune-expressed\n(HPA, nTPM≥1)',
              'Druggable genome\n(Finan et al. 2017)']

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.6, 5.0),
                                   gridspec_kw={'width_ratios': [1.05, 1]})

    # ---------- Panel A: rate bars ----------
    iei_rate = [sub.loc[k, 'p_tgt_iei'] * 100 for k in order3]
    non_rate = [sub.loc[k, 'p_tgt_non'] * 100 for k in order3]
    x = np.arange(len(order3)); w = 0.36
    axA.bar(x - w/2, iei_rate, w, label='IEI-associated genes', color='#c0392b')
    axA.bar(x + w/2, non_rate, w, label='all other genes', color='#9aa5a8')
    axA.set_xticks(x); axA.set_xticklabels(disp3, fontsize=7.5)
    axA.set_ylabel('% that are targets of an\nFDA-approved immune drug', labelpad=10)
    axA.set_ylim(0, 40)
    set_frame(axA)
    axA.legend(frameon=False, loc='upper left', fontsize=8, bbox_to_anchor=(0.0, 1.0))
    for xi, v in zip(x - w/2, iei_rate):
        axA.annotate(f'{v:.1f}', (xi, v + 0.7), ha='center', fontsize=8,
                     fontweight='bold', color='#c0392b')
    for xi, v in zip(x + w/2, non_rate):
        axA.annotate(f'{v:.1f}', (xi, v + 0.7), ha='center', fontsize=8, color='#555')
    # fold-enrichment bracket above each pair — full-precision ratio, 1 decimal
    for xi, i, n in zip(x, iei_rate, non_rate):
        top = max(i, n) + 4.0
        axA.annotate('', xy=(xi - w/2, top), xytext=(xi + w/2, top),
                     arrowprops=dict(arrowstyle='-', color='#2c3e50', lw=0.9))
        axA.annotate(f'{i / n:.1f}×', (xi, top + 0.3), ha='center', fontsize=8,
                     color='#2c3e50', fontweight='bold')
    panel_letter(axA, 'a')

    # ---------- Panel B: odds-ratio forest ----------
    ors = np.array([sub.loc[k, 'OR']    for k in order5])
    los = np.array([sub.loc[k, 'OR_lo'] for k in order5])
    his = np.array([sub.loc[k, 'OR_hi'] for k in order5])
    y = np.arange(len(order5))[::-1]
    axB.errorbar(ors, y, xerr=[ors - los, his - ors], fmt='o', color='#2c3e50',
                 ms=8, capsize=3.5, lw=1.6)
    axB.axvline(1, color='#c0392b', ls='--', lw=1.2)
    axB.text(1.08, y.max() + 0.35, 'no enrichment', color='#c0392b', fontsize=7.5, va='center')
    for yi, o, lo, hi in zip(y, ors, los, his):
        axB.annotate(f'{o:.1f} ({lo:.1f}–{hi:.1f})', (hi + 0.2, yi), va='center', fontsize=7.8)
    axB.set_yticks(y); axB.set_yticklabels(disp5, fontsize=7.8)
    axB.set_xlabel('Odds ratio (95% CI)\nenrichment of drug-target status among IEI genes', labelpad=6)
    axB.set_xlim(0, 9.4); axB.set_ylim(-0.75, len(order5) - 0.15)
    set_frame(axB)
    panel_letter(axB, 'b')

    fig.suptitle('IEI-associated genes are enriched for FDA immune-drug targets '
                 'across every background universe',
                 fontsize=10, y=1.01, x=0.02, ha='left')
    fig.text(0.5, -0.06,
             'Immune-expressed universe: Human Protein Atlas, 19 sorted immune cell types '
             '(max nTPM ≥ threshold). Druggable genome: Finan et al. 2017. Enrichment by '
             'two-sided Fisher exact test; all p < 1×10⁻¹⁰. n = 505 IEI genes, '
             '723 FDA immune-drug target genes, 72 shared.',
             ha='center', fontsize=6.6, color='#555', wrap=True)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches='tight')
    return fig


if __name__ == "__main__":
    import sys
    from pathlib import Path
    _REPO = Path(__file__).resolve().parents[2]
    # build_enrichment_table lives in the sibling stage folder src/iei_enrichment/.
    sys.path.insert(0, str(_REPO / "src" / "iei_enrichment"))
    from make_enrichment_table import build_enrichment_table
    res, _ = build_enrichment_table()
    out = _REPO / "plots" / "iei_enrichment_manuscript.png"
    make_figure(res, outfile=str(out))
    print(f"wrote {out}")
