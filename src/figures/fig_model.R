#!/usr/bin/env Rscript
# make_model_fig.R
#
# Two-panel integrative PU-model figure (ggplot2 / patchwork) reproducing the
# figure that make_model_fig.py draws. This script reads ONLY the model run's
# already-written outputs (nominations.tsv, metrics.json) -- it does NOT re-run
# the model. Theme mirrors Part 2's make_regulator_burden_figure.R.
#
#   Panel a  PU score distributions -- known immune targets (label_immune==1) vs
#                                       unlabeled genes, log-count y-axis, with the
#                                       three highlighted candidates (MAP3K14,
#                                       RORC, MALT1) called out.
#   Panel b  GWAS ablation lollipop -- grouped-CV enrichment@100 with and without
#                                       the Open Targets GWAS score, showing
#                                       genetics carries about half the signal.
#
# Inputs (path resolved relative to this script; the model outputs live in the
# repo-level outputs/ dir, one level up from Part 5):
#   ../outputs/nominations.tsv   columns: gene, pu_score, label_immune, ...
#   ../outputs/metrics.json      ablations[].variant + ablations[]."enrichment@100"
#
# Output (written next to this script):
#   model_manuscript.png   (300 dpi; does not overwrite the Python-made
#                           model_fig.png)
#
# Usage:  Rscript make_model_fig.R

# Force a UTF-8 locale so Unicode glyphs in labels (× − ↑) render instead of
# being mangled to one dot per byte under a C locale.
for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

# --- Resolve paths relative to this script's location -----------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()

need <- function(f) {
  p <- normalizePath(file.path(script_dir, f), mustWork = FALSE)
  if (!file.exists(p))
    stop("Missing ", p,
         "\nRun `python -m model.run` first to generate the model outputs.")
  p
}

# --- Shared style (mirrors Part 2's make_regulator_burden_figure.R) ---------
CANDIDATES <- c("MAP3K14", "RORC", "MALT1")
col_unl    <- "#b0b0b0"   # unlabeled genes (META_GREY)
col_known  <- "#2a7ab0"   # known immune targets / with-GWAS (blue)
col_ctrl   <- "#c0c0c0"   # without-GWAS (grey)
col_cand   <- "#d1495b"   # candidate call-outs (red)

# Font with full Unicode coverage so × − ↑ render (macOS "sans" drops them).
font_family <- "Arial Unicode MS"
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

base_theme <- theme_classic(base_size = 9) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = 9),
        axis.text   = element_text(size = 7),
        plot.title  = element_text(size = 8.5, hjust = 0),
        plot.margin = margin(6, 10, 6, 10))

# --- Load the model outputs -------------------------------------------------
nom <- read.delim(need(file.path("..", "..", "outputs", "nominations.tsv")),
                  stringsAsFactors = FALSE)
metrics <- jsonlite::fromJSON(need(file.path("..", "..", "outputs", "metrics.json")),
                              simplifyDataFrame = TRUE)

# ===================== Panel a: PU score distributions =====================
# Precompute histogram counts per group (hist plot=FALSE) so empty bins are
# simply dropped rather than becoming log(0) = -Inf on the log y-axis.
brks <- seq(0, 1, length.out = 46)
hist_df <- function(x, grp) {
  h <- hist(x, breaks = brks, plot = FALSE)
  data.frame(mid = h$mids, count = h$counts, group = grp)
}
pos <- nom$pu_score[nom$label_immune == 1]
unl <- nom$pu_score[nom$label_immune == 0]
bars <- rbind(
  hist_df(unl, sprintf("Unlabeled genes (n=%s)",      format(length(unl), big.mark = ","))),
  hist_df(pos, sprintf("Known immune targets (n=%s)", format(length(pos), big.mark = ",")))
)
bars <- bars[bars$count > 0, ]
grp_levels <- c(sprintf("Unlabeled genes (n=%s)",      format(length(unl), big.mark = ",")),
                sprintf("Known immune targets (n=%s)", format(length(pos), big.mark = ",")))
bars$group <- factor(bars$group, levels = grp_levels)

# candidate scores, highest at top for the stacked call-out labels
cand <- data.frame(gene = CANDIDATES,
                   pu_score = vapply(CANDIDATES,
                                     function(g) nom$pu_score[nom$gene == g][1], numeric(1)))
cand <- cand[order(-cand$pu_score), ]
cand$ylab <- 10 ^ (3.55 - 0.62 * (seq_len(nrow(cand)) - 1))   # stagger label heights

panelA <- ggplot(bars, aes(mid, count, fill = group)) +
  geom_col(width = diff(brks)[1], position = "identity", alpha = 0.8) +
  # vertical mark at each candidate's score (mirrors axvline(s, ymax=0.52))
  geom_segment(data = cand, inherit.aes = FALSE,
               aes(x = pu_score, xend = pu_score, y = 0.7, yend = 10^2.3),
               colour = col_cand, linewidth = 0.4) +
  # thin connector from just RIGHT of the offset label to the mark (mirrors the
  # .py's annotate arrowprops; start past the text so it doesn't cover the label)
  geom_segment(data = cand, inherit.aes = FALSE,
               aes(x = 0.47, xend = pu_score, y = ylab, yend = 10^2.3),
               colour = col_cand, linewidth = 0.25) +
  geom_text(data = cand, inherit.aes = FALSE,
            aes(x = 0.31, y = ylab, label = gene),
            colour = col_cand, fontface = "italic", size = 2.4, hjust = 0) +
  scale_fill_manual(values = setNames(c(col_unl, col_known), grp_levels), name = NULL) +
  scale_y_log10(limits = c(0.7, 3e4),
                breaks = c(1, 10, 100, 1000, 10000),
                labels = c("1", "10", "100", "1000", "10000")) +
  labs(x = "PU nomination score", y = "Number of genes",
       title = "Model separates known targets from the genome", tag = "a") +
  base_theme +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.5, 0.98),
        legend.justification   = c(0.5, 1),
        legend.text            = element_text(size = 6.8, family = font_family),
        legend.key.size        = unit(9, "pt"),
        plot.tag               = element_text(size = 13))

# ==================== Panel b: GWAS ablation (enrichment@100) ====================
abl     <- metrics$ablations
enrich  <- abl[["enrichment@100"]]
get_enr <- function(v) enrich[abl$variant == v][1]
# numeric x throughout so the numeric annotate() call and the lollipop share one
# continuous x scale (mixing factor geoms with numeric annotate() errors out).
gw <- data.frame(
  x      = c(1, 2),
  label  = c("+ GWAS score", "– GWAS score"),
  enrich = c(get_enr("(i)  with gwas_score"), get_enr("(i)  without gwas_score")),
  fill   = c(col_known, col_ctrl)
)

panelB <- ggplot(gw, aes(x, enrich)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +                                   # no-enrichment ref
  annotate("text", x = 2.5, y = 1, label = "no enrichment", vjust = -0.5,
           hjust = 1, size = 2.3, colour = "grey45") +
  geom_segment(aes(xend = x, y = 0, yend = enrich, colour = fill),
               linewidth = 0.7) +                                  # lollipop stem
  geom_point(aes(colour = fill), size = 4.2) +
  geom_text(aes(label = sprintf("%.1f×", enrich)), vjust = -1.1, size = 3) +
  annotate("text", x = 2.5, y = 9.6, label = "↑ higher = better",
           hjust = 1, size = 2.3, colour = "grey50") +
  scale_colour_identity() +
  scale_x_continuous(breaks = gw$x, labels = gw$label, limits = c(0.5, 2.6)) +
  scale_y_continuous(limits = c(0, 9.8), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Enrichment in top 100 (grouped-CV)",
       title = "Genetics carries about half the signal", tag = "b") +
  base_theme +
  theme(plot.tag = element_text(size = 13))

# ============================== Compose + save ==============================
fig <- panelA + panelB + plot_layout(widths = c(1.35, 1))

outfile <- file.path(script_dir, "..", "..", "plots", "model_manuscript.png")
# ragg's AGG device renders the Unicode glyphs the base PNG device drops;
# fall back to cairo if ragg is unavailable.
dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "cairo"
ggsave(outfile, fig, width = 9.4, height = 4.3, dpi = 300, bg = "white", device = dev)
cat("wrote", outfile, "\n")
