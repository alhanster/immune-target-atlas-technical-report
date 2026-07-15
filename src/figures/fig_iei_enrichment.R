#!/usr/bin/env Rscript
# make_manuscript_figure.R
#
# Two-panel manuscript figure (ggplot2 / patchwork port of make_manuscript_figure.py):
# enrichment of IEI-associated genes among approved immune-drug targets, across five
# background universes.
#
#   Panel A  rate bars       -- % of IEI vs non-IEI genes that are targets, in three
#                               headline universes, with fold-enrichment brackets
#   Panel B  odds-ratio forest -- OR + 95% CI across all five universes
#
# Input (path resolved relative to this script's location):
#   iei_enrichment_by_universe.csv
#       Written by build_enrichment_table() in make_enrichment_table.py. One row per
#       background universe with columns:
#         universe, N, a, b, c, d, p_tgt_iei, p_tgt_non, RR, OR, OR_lo, OR_hi, p
#
# Output (written next to this script):
#   iei_enrichment_manuscript.png   (300 dpi)
#
# Usage:  Rscript make_manuscript_figure.R

# Force a UTF-8 locale so the Unicode glyphs in labels (≥ × ∩ – ⁻¹⁰) render
# instead of being mangled to one dot per byte under a C locale.
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

csv_path <- file.path(script_dir, "..", "..", "data", "derived", "iei_enrichment_by_universe.csv")
if (!file.exists(csv_path)) {
  stop("Missing ", csv_path,
       "\nRun `python3 make_enrichment_table.py` first to generate it.")
}
res <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(res) <- res$universe

# --- Universe orderings + display labels (mirror the Python script) ---------
order5 <- c("All genes (HPA-20k proxy)", "Immune-expressed nTPM>=1",
            "Immune-expressed nTPM>=5", "Druggable genome",
            "Druggable ∩ immune-expr(>=1)")
disp5  <- c("All genes (~20,000)", "Immune-expressed (HPA, nTPM≥1)",
            "Immune-expressed (HPA, nTPM≥5)", "Druggable genome (Finan 2017)",
            "Druggable ∩ immune-expressed")
order3 <- c("All genes (HPA-20k proxy)", "Immune-expressed nTPM>=1", "Druggable genome")
disp3  <- c("All genes\n(~20,000)", "Immune-expressed\n(HPA, nTPM≥1)",
            "Druggable genome\n(Finan et al. 2017)")

col_iei   <- "#56B4E9"   # IEI-associated bars
col_other <- "#999999"   # non-IEI bars
col_dark  <- "#2c3e50"   # forest points / lines

# Fisher-exact significance stars from a p-value.
sig <- function(p) ifelse(p < .001, "***",
                   ifelse(p < .01,  "**",
                   ifelse(p < .05,  "*", "ns")))

# Use a font with full Unicode coverage so ≥ × ∩ – ⁻¹⁰ render (macOS default
# "sans" drops them). Applied to theme text AND geom_text/label defaults.
font_family <- "Arial Unicode MS"
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

base_theme <- theme_classic(base_size = 9) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = 9),
        axis.text   = element_text(size = 7),
        plot.margin = margin(6, 10, 6, 10))

# ============================ Panel A: rate bars ============================
iei_rate <- res[order3, "p_tgt_iei"] * 100
non_rate <- res[order3, "p_tgt_non"] * 100
xn <- seq_along(order3)
w  <- 0.36

barsA <- data.frame(
  x     = c(xn - w / 2, xn + w / 2),
  rate  = c(iei_rate, non_rate),
  group = factor(rep(c("IEI-associated genes", "all other genes"), each = length(xn)),
                 levels = c("IEI-associated genes", "all other genes"))
)
# fold-enrichment bracket geometry (one per universe pair)
brk <- data.frame(
  xl   = xn - w / 2, xr = xn + w / 2, xc = xn,
  top  = pmax(iei_rate, non_rate) + 4.0,
  fold = sprintf("%.1f×", iei_rate / non_rate),
  star = sig(res[order3, "p"])
)

panelA <- ggplot(barsA, aes(x, rate, fill = group)) +
  geom_col(width = w, position = position_identity(), alpha = 0.7) +
  scale_fill_manual(values = c("IEI-associated genes" = col_iei,
                               "all other genes"       = col_other), name = NULL) +
  # per-bar value labels (black, both series)
  geom_text(data = data.frame(x = xn - w / 2, rate = iei_rate),
            aes(x, rate + 0.7, label = sprintf("%.1f", rate)),
            inherit.aes = FALSE, size = 2.9, colour = "black") +
  geom_text(data = data.frame(x = xn + w / 2, rate = non_rate),
            aes(x, rate + 0.7, label = sprintf("%.1f", rate)),
            inherit.aes = FALSE, size = 2.9, colour = "black") +
  # fold-enrichment brackets, significance stars, then fold labels (black)
  geom_segment(data = brk, aes(x = xl, xend = xr, y = top, yend = top),
               inherit.aes = FALSE, colour = col_dark, linewidth = 0.35) +
  geom_text(data = brk, aes(x = xc, y = top + 1.0, label = star),
            inherit.aes = FALSE, size = 3.4, colour = "black") +
  geom_text(data = brk, aes(x = xc, y = top + 2.6, label = fold),
            inherit.aes = FALSE, size = 2.9, colour = "black") +
  scale_x_continuous(breaks = xn, labels = disp3) +
  scale_y_continuous(limits = c(0, 44), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL,
       y = "% that are targets of an\napproved immune drug",
       tag = "a") +
  base_theme +
  theme(axis.text.x            = element_text(size = 7.8),
        legend.position        = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification   = c(0, 1),
        legend.text            = element_text(size = 8),
        legend.key.size        = unit(9, "pt"),
        plot.tag               = element_text(size = 13))

# ========================= Panel B: odds-ratio forest =========================
forest <- data.frame(
  y     = rev(seq_along(order5)),   # first universe at the top (y = 5)
  OR    = res[order5, "OR"],
  lo    = res[order5, "OR_lo"],
  hi    = res[order5, "OR_hi"]
)
forest$txt <- sprintf("%.1f (%.1f–%.1f) %s",
                      forest$OR, forest$lo, forest$hi, sig(res[order5, "p"]))

panelB <- ggplot(forest, aes(OR, y)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#c0392b", linewidth = 0.5) +
  annotate("text", x = 1.08, y = max(forest$y) + 0.35, label = "No Enrichment",
           colour = "#c0392b", size = 2.9, hjust = 0) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = col_dark,
                 linewidth = 0.6) +
  geom_point(colour = col_dark, size = 2.6) +
  geom_text(aes(x = hi + 0.2, label = txt), hjust = 0, size = 2.9) +
  scale_x_continuous(limits = c(0, 10.2), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = forest$y, labels = disp5,
                     limits = c(0.25, length(order5) + 0.85)) +
  coord_cartesian(clip = "off") +
  labs(x = "Odds ratio (95% CI)\nenrichment of drug-target status among IEI genes",
       y = NULL, tag = "b") +
  base_theme +
  theme(axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(size = 7.8),
        plot.tag     = element_text(size = 13))

# ============================== Compose + save ==============================
fig <- panelA + panelB +
  plot_layout(widths = c(1.05, 1))

outfile <- file.path(script_dir, "..", "..", "plots", "iei_enrichment_manuscript.png")
# ragg's AGG device renders the Unicode glyphs (≥ × ∩ – ⁻¹⁰) that the base PNG
# device drops; fall back to cairo if ragg is unavailable.
dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "cairo"
ggsave(outfile, fig, width = 11.6, height = 5.5, dpi = 300, bg = "white", device = dev)
cat("wrote", outfile, "\n")
