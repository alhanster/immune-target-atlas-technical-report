#!/usr/bin/env Rscript
# make_regulator_burden_figure.R
#
# Two-panel manuscript figure (ggplot2 / patchwork) reproducing the figure that
# regulator_burden_pipeline.py draws in make_figure(). This script reads ONLY
# the pipeline's already-written output CSVs -- it does NOT re-run the pipeline
# (no S3 stream, no h5py). Theme mirrors Part 1's make_manuscript_figure.R.
#
#   Panel a  signed-QQ           -- observed vs expected signed -log10 P of the
#                                   genome-wide regulator-burden correlation,
#                                   across the three culture conditions
#   Panel b  IEI core-gene bars  -- top-18 IEI genes nominated as core genes for
#                                   lymphocyte count (-log10 P, colored by condition)
#
# Inputs (path resolved relative to this script; must already exist):
#   regulator_burden_scores_Rest.csv     } column `signed_log10p`  (Panel a)
#   regulator_burden_scores_Stim8hr.csv  }
#   regulator_burden_scores_Stim48hr.csv }
#   regulator_burden_scores_all.csv        columns gene, condition, p_beta,
#                                          signed_log10p              (Panel b)
#   ../../Data/Gene List/IEI_gene_list.csv    IEI panel; Panel b keeps only these
#                                          genes ('Gene' column)
#
# Output (written next to this script):
#   regulator_burden_manuscript.png   (300 dpi; does not overwrite the
#                                      Python-made regulator_burden_figure.png)
#
# Usage:  Rscript make_regulator_burden_figure.R

# Force a UTF-8 locale so Unicode glyphs in labels (− × ≥) render instead of
# being mangled to one dot per byte under a C locale.
for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# --- Resolve paths relative to this script's location -----------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()

need <- function(f) {
  p <- file.path(script_dir, "..", "..", "data", "derived", f)
  if (!file.exists(p))
    stop("Missing ", p,
         "\nRun `python regulator_burden_pipeline.py` first to generate outputs.")
  p
}

# --- Shared style (mirrors Part 1's make_manuscript_figure.R) ---------------
conds    <- c("Rest", "Stim8hr", "Stim48hr")
cond_col <- c(Rest = "#4393C3", Stim8hr = "#D6604D", Stim48hr = "#E08214")  # Python colmap

# Font with full Unicode coverage so − × ≥ render (macOS "sans" drops them).
font_family <- "Arial Unicode MS"
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

base_theme <- theme_classic(base_size = 9) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = 9),
        axis.text   = element_text(size = 7),
        plot.title  = element_text(size = 8.5, hjust = 0.5),
        plot.margin = margin(6, 10, 6, 10))

# --- Genome-wide per-gene scores (read once; used by Panel a) ---------------
scores_long <- do.call(rbind, lapply(conds, function(C) {
  d <- read.csv(need(sprintf("regulator_burden_scores_%s.csv", C)))
  data.frame(condition = C, gene = d$gene, signed_log10p = d$signed_log10p)
}))
scores_long$condition <- factor(scores_long$condition, levels = conds)

# ===================== Panel a: signed-QQ across conditions =====================
# For each condition, contrast observed signed -log10 P against the expectation
# under a uniform null (same transform as make_figure() in the .py).
qq <- do.call(rbind, lapply(conds, function(C) {
  obs <- sort(scores_long$signed_log10p[scores_long$condition == C])
  n   <- length(obs)
  q   <- (seq_len(n) - 0.5) / n
  exp <- sign(q - 0.5) * -log10(1 - 2 * abs(q - 0.5) + 1e-300)
  data.frame(condition = C, exp = exp, obs = obs)
}))
qq$condition <- factor(qq$condition, levels = conds)

panelA <- ggplot(qq, aes(exp, obs, colour = condition)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40", linewidth = 0.35) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_colour_manual(values = cond_col, name = NULL) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-8, 8)) +
  labs(x = expression("Expected signed " * -log[10] * P),
       y = expression("Observed signed " * -log[10] * P),
       title = "Regulator-burden correlation vs lymphocyte count",
       tag = "a") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  base_theme +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification   = c(0, 1),
        legend.text            = element_text(size = 7, family = font_family),
        legend.key.size        = unit(9, "pt"),
        plot.tag               = element_text(size = 13))

# ==================== Panel b: IEI genes nominated as core genes ====================
# IEI selection lives here (not in the pipeline): read the IEI panel and keep
# only those genes from the genome-wide scores. Panel is the Part 1 copy at
# ../../Data/Gene List/IEI_gene_list.csv (identical to the vendored inputs/ copy).
iei_csv <- file.path(script_dir, "..", "..", "data", "reference", "IEI_gene_list.csv")
if (!file.exists(iei_csv)) stop("IEI gene list not found: ", iei_csv)
iei_genes <- read.csv(iei_csv, stringsAsFactors = FALSE)$Gene
iei <- read.csv(need("regulator_burden_scores_all.csv"), stringsAsFactors = FALSE)
iei <- iei[iei$gene %in% iei_genes, ]              # keep IEI genes only
iei$neglogp <- -log10(iei$p_beta)                  # per (gene, condition)

# Rank IEI genes by the largest absolute Rest-vs-stim gap in -log10P (Rest vs
# Stim8hr, or Rest vs Stim48hr, whichever is larger); keep the top 15, then show
# every condition's value as a colored dot.
wide  <- reshape(iei[, c("gene", "condition", "neglogp")],
                 idvar = "gene", timevar = "condition", direction = "wide")
gdiff <- pmax(abs(wide$neglogp.Rest - wide$neglogp.Stim8hr),
              abs(wide$neglogp.Rest - wide$neglogp.Stim48hr), na.rm = TRUE)
names(gdiff) <- wide$gene
top_genes  <- names(sort(gdiff, decreasing = TRUE))[1:15]
dot        <- iei[iei$gene %in% top_genes, ]
dot$gene   <- factor(dot$gene, levels = rev(top_genes))       # largest diff on top
dot$condition <- factor(dot$condition, levels = conds)

panelB <- ggplot(dot, aes(x = neglogp, y = gene)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             colour = "grey55", linewidth = 0.3) +            # P = 0.05 reference
  geom_line(aes(group = gene), colour = "grey80", linewidth = 0.4) +  # connect dots
  geom_point(aes(colour = condition), size = 2.2) +
  scale_colour_manual(values = cond_col, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = expression(-log[10] * P * "  (regulator-burden)"),
       y = NULL,
       title = "IEI genes nominated as core genes\nfor lymphocyte count",
       tag = "b") +
  base_theme +
  theme(axis.text.y            = element_text(size = 7.5, face = "italic"),
        legend.position        = "inside",
        legend.position.inside = c(0.98, 0.02),
        legend.justification   = c(1, 0),
        legend.text            = element_text(size = 7, family = font_family),
        legend.key.size        = unit(9, "pt"),
        plot.tag               = element_text(size = 13))

# ============================== Compose + save ==============================
fig <- panelA + panelB + plot_layout(widths = c(1, 1))

outfile <- file.path(script_dir, "..", "..", "plots", "regulator_burden_manuscript.png")
# ragg's AGG device renders the Unicode glyphs the base PNG device drops;
# fall back to cairo if ragg is unavailable.
dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "cairo"
ggsave(outfile, fig, width = 9.6, height = 4.4, dpi = 300, bg = "white", device = dev)
cat("wrote", outfile, "\n")
