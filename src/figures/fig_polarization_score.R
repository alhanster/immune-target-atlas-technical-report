#!/usr/bin/env Rscript
# make_polarization_score_figure.R
#
# Two-panel diagnostic figure (ggplot2 / patchwork) reproducing the figure that
# polarization_score.py draws in make_figure(). This script reads ONLY the
# already-written output CSV (polarization_score.csv) -- it does NOT re-run the
# analysis (no coefficient CSV, no scipy, no score recomputation). Theme mirrors
# make_regulator_burden_figure.R.
#
#   Panel a  8h vs 48h z scatter -- reinforcement on the diagonal, cancellation
#                                   off it; known Th1/Th2 regulators highlighted
#   Panel b  signed-score hist   -- Th2 (negative) vs Th1 (positive) poles
#
# Input (path resolved relative to this script; must already exist):
#   polarization_score.csv   columns gene, polarization_score, signs_agree,
#                            z_Stim8hr, z_Stim48hr, known_regulator,
#                            regulator_type
#
# Output (written next to this script):
#   polarization_score_manuscript.png   (300 dpi; does not overwrite the
#                                        Python-made polarization_score_diagnostics.png)
#
# Usage:  Rscript make_polarization_score_figure.R

# Force a UTF-8 locale so Unicode glyphs in labels (− √ ₈) render instead of
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
         "\nRun `python polarization_score.py` first to generate outputs.")
  p
}

# --- Shared style (mirrors make_regulator_burden_figure.R) ------------------
# Palette carried over verbatim from make_figure() in polarization_score.py.
META_GREY  <- "#888888"                        # neutral reference-line grey
col_cancel <- "#d9b38c"; col_reinf <- "#2a7ab0"
col_th2    <- "#c1666b"; col_th1   <- "#3b5b78"
col_hist   <- "#c9ccd1"
cat_lvls   <- c("signs disagree (cancel)", "signs agree (reinforce)",
                "Known Th2 regulator", "Known Th1 regulator")
cat_cols   <- setNames(c(col_cancel, col_reinf, col_th2, col_th1), cat_lvls)

# Font with full Unicode coverage so − √ render (macOS "sans" drops them).
font_family <- "Arial Unicode MS"
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

base_theme <- theme_classic(base_size = 9) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = 9),
        axis.text   = element_text(size = 7),
        plot.title  = element_text(size = 8.5, hjust = 0.5),
        plot.margin = margin(6, 10, 6, 10))

# --- Per-gene scores (read once; used by both panels) -----------------------
S <- read.csv(need("polarization_score.csv"), stringsAsFactors = FALSE,
              row.names = NULL)
# pandas writes booleans as "True"/"False" -- R's as.logical() maps those to NA,
# so normalise explicitly.
as_bool <- function(x) toupper(as.character(x)) == "TRUE"
S$signs_agree     <- as_bool(S$signs_agree)
S$known_regulator <- as_bool(S$known_regulator)

kn  <- S[S$known_regulator, ]
th2 <- kn[kn$regulator_type == "Th2 regulator" & !is.na(kn$regulator_type), ]
th1 <- kn[kn$regulator_type == "Th1 regulator" & !is.na(kn$regulator_type), ]

# =============== Panel a: 8h vs 48h z (Stouffer reinforcement) ===============
# Background points split by sign agreement, drawn first; known Th1/Th2
# regulators overplotted larger on top. Constant colour aes per layer so the
# four categories accumulate into a single legend (as in make_figure()).
dis <- S[!S$signs_agree, ]
ag  <- S[ S$signs_agree, ]

L <- 75
panelA <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, colour = META_GREY) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = META_GREY) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = META_GREY) +
  geom_point(data = dis, aes(z_Stim8hr, z_Stim48hr, colour = cat_lvls[1]),
             size = 0.45, alpha = 0.50) +
  geom_point(data = ag,  aes(z_Stim8hr, z_Stim48hr, colour = cat_lvls[2]),
             size = 0.45, alpha = 0.50) +
  geom_point(data = th2, aes(z_Stim8hr, z_Stim48hr, colour = cat_lvls[3]),
             size = 2.4) +
  geom_point(data = th1, aes(z_Stim8hr, z_Stim48hr, colour = cat_lvls[4]),
             size = 2.4) +
  scale_colour_manual(values = cat_cols, breaks = cat_lvls, name = NULL) +
  coord_cartesian(xlim = c(-L, L), ylim = c(-L, L)) +
  labs(x = expression("polarization z at Stim8hr"),
       y = expression("polarization z at Stim48hr"),
       title = expression("Score = (z"[8*hr] * " + z"[48*hr] * ")/" * sqrt(2) *
                          "  — Stouffer combination"),
       tag = "a") +
  guides(colour = guide_legend(override.aes =
           list(size = c(1.3, 1.3, 2.4, 2.4), alpha = 1))) +
  base_theme +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.98, 0.02),
        legend.justification   = c(1, 0),
        legend.text            = element_text(size = 7, family = font_family),
        legend.key.size        = unit(9, "pt"),
        plot.tag               = element_text(size = 13))

# Anchor-gene labels: nudges are the .py's point offsets converted to data
# units (~0.5 unit/pt at this scale), keeping the same placement.
lab_nudge <- list(IL4R  = c( 3, -1), STAT6  = c(  3, -1), GATA3 = c(3,  1.5),
                  TBX21 = c( 3,  2), IFNGR1 = c(-21,  1), RARA  = c(3, -4))
lab_df <- do.call(rbind, lapply(names(lab_nudge), function(g) {
  r <- S[S$gene == g, ]
  if (nrow(r) == 0) return(NULL)
  data.frame(gene = g, x = r$z_Stim8hr[1], y = r$z_Stim48hr[1],
             nx = lab_nudge[[g]][1], ny = lab_nudge[[g]][2])
}))
if (!is.null(lab_df))
  panelA <- panelA +
    geom_text(data = lab_df, aes(x + nx, y + ny, label = gene),
              size = 2.1, hjust = 0, colour = "grey20")

# ==================== Panel b: signed-score distribution ====================
brks <- seq(min(S$polarization_score), max(S$polarization_score), length.out = 81)
panelB <- ggplot(S, aes(polarization_score)) +
  geom_histogram(breaks = brks, fill = col_hist, alpha = 0.8) +
  geom_vline(data = th2, aes(xintercept = polarization_score),
             colour = col_th2, linewidth = 0.3, alpha = 0.8) +
  geom_vline(data = th1, aes(xintercept = polarization_score),
             colour = col_th1, linewidth = 0.3, alpha = 0.8) +
  geom_vline(xintercept = 0, colour = "#222222", linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(x = expression("polarization score  (signed Stouffer z)"),
       y = "number of genes",
       title = "Signed score separates poles:\nKnown Th2 (red) negative, Known Th1 (blue) positive",
       tag = "b") +
  base_theme +
  theme(plot.tag = element_text(size = 13))

# ============================== Compose + save ==============================
fig <- panelA + panelB + plot_layout(widths = c(1, 1))

outfile <- file.path(script_dir, "..", "..", "plots", "polarization_score_manuscript.png")
# ragg's AGG device renders the Unicode glyphs the base PNG device drops;
# fall back to cairo if ragg is unavailable.
dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "cairo"
ggsave(outfile, fig, width = 10.5, height = 4.6, dpi = 300, bg = "white", device = dev)
cat("wrote", outfile, "\n")
