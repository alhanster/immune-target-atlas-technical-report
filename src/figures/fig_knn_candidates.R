#!/usr/bin/env Rscript
# make_knn_candidates_figure.R
#
# Signature-embedding manuscript figure (ggplot2 / patchwork). Reads the
# per-condition signature-PC tables written by embed_signatures.py, computes a
# UMAP per condition (uwot), and draws a 2x2 grid:
#   rows    = culture condition (Stim8hr, Stim48hr)
#   columns = coloring: FDA membership | nearest-target hubs
# to show whether FDA drug-target anchors form coherent mechanistic islands in
# knockdown-signature space. Theme mirrors Part 2's make_regulator_burden_figure.R.
#
# Inputs (path resolved relative to this script; run embed_signatures.py first):
#   knn_signature_pcs_Stim8hr.csv   } columns: gene, is_fda_target,
#   knn_signature_pcs_Stim48hr.csv  } nearest_fda_target, PC1..PC50
#
# Output (written next to this script):
#   knn_immune_target_manuscript.png   (300 dpi)
#
# Usage:  Rscript make_knn_candidates_figure.R

# Force a UTF-8 locale so Unicode glyphs in labels (· −) render instead of
# being mangled to one dot per byte under a C locale.
for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(uwot)
})

# --- Resolve paths relative to this script's location -----------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()

need <- function(f) {
  p <- file.path(script_dir, "..", "..", "data", "derived", f)
  if (!file.exists(p))
    stop("Missing ", p,
         "\nRun `python embed_signatures.py` first to generate the PC tables.")
  p
}

# --- Shared style (mirrors Part 2's make_regulator_burden_figure.R) ---------
conds  <- c("Stim8hr", "Stim48hr")
N_PCS  <- 50
N_HUBS <- 8

col_cand   <- "#BDC3C7"   # non-FDA candidate / "other" (gray)
col_anchor <- "#2E86C1"   # FDA anchor (blue)
# Okabe-Ito colorblind-safe palette for the 8 nearest-target hubs.
okabe_ito  <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Font with full Unicode coverage so · − render (macOS "sans" drops them).
font_family <- "Arial Unicode MS"
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

base_theme <- theme_classic(base_size = 9) +
  theme(text         = element_text(family = font_family),
        axis.title   = element_text(size = 8),
        axis.text    = element_text(size = 6),
        plot.title   = element_text(size = 8.5, hjust = 0.5),
        plot.margin  = margin(6, 10, 6, 10),
        legend.text  = element_text(size = 6,   family = font_family),
        legend.title = element_text(size = 6.5, family = font_family),
        legend.key.size = unit(8, "pt"))

# --- Per-condition UMAP embedding of the signature PCs ----------------------
embed <- function(C) {
  d  <- read.csv(need(sprintf("knn_signature_pcs_%s.csv", C)), stringsAsFactors = FALSE)
  pc <- as.matrix(d[, sprintf("PC%d", seq_len(N_PCS))])
  set.seed(1)                                   # reproducible UMAP (+ n_threads=1)
  um <- uwot::umap(pc, n_neighbors = 15, min_dist = 0.1,
                   metric = "euclidean", n_threads = 1)
  d$UMAP1 <- um[, 1]; d$UMAP2 <- um[, 2]
  d$is_fda <- as.logical(d$is_fda_target)
  d
}
emb <- setNames(lapply(conds, embed), conds)

# --- Panel builders ---------------------------------------------------------
# Left column: FDA anchors (blue) over non-FDA candidates (gray) -- do the
# known drug targets occupy coherent regions, or scatter diffusely?
panel_membership <- function(d, C) {
  ggplot(d, aes(UMAP1, UMAP2)) +
    geom_point(data = d[!d$is_fda, ], colour = col_cand, size = 0.5, alpha = 0.5) +
    geom_point(data = d[d$is_fda, ],  fill = col_anchor, shape = 21,
               size = 1.5, colour = "black", stroke = 0.15) +
    labs(x = "UMAP1", y = "UMAP2", title = sprintf("%s · FDA membership", C)) +
    base_theme
}

# Right column: color every gene by which of the 8 most-recurrent FDA targets it
# is nearest to (its mechanistic "island"); rest gray. Labels at island medians.
panel_hubs <- function(d, C) {
  hubs <- names(sort(table(d$nearest_fda_target), decreasing = TRUE))[seq_len(N_HUBS)]
  d$hub <- factor(ifelse(d$nearest_fda_target %in% hubs, d$nearest_fda_target, "other"),
                  levels = c(hubs, "other"))
  pal  <- setNames(c(okabe_ito, col_cand), c(hubs, "other"))
  cent <- do.call(rbind, lapply(hubs, function(h) {
    s <- d[d$hub == h, ]
    data.frame(hub = h, UMAP1 = median(s$UMAP1), UMAP2 = median(s$UMAP2))
  }))
  ggplot(d, aes(UMAP1, UMAP2)) +
    geom_point(data = d[d$hub == "other", ], colour = col_cand, size = 0.5, alpha = 0.4) +
    geom_point(data = d[d$hub != "other", ], aes(fill = hub), shape = 21,
               size = 1.5, colour = "black", stroke = 0.15) +
    geom_text(data = cent, aes(label = hub), fontface = "italic", size = 2, colour = "black") +
    scale_fill_manual(values = pal, breaks = hubs, name = "nearest FDA target") +
    labs(x = "UMAP1", y = "UMAP2", title = sprintf("%s · nearest-target hubs", C)) +
    base_theme
}

# --- Compose 2x2 (rows = condition, cols = coloring) ------------------------
fig <-
  (panel_membership(emb[["Stim8hr"]],  "Stim8hr")  | panel_hubs(emb[["Stim8hr"]],  "Stim8hr")) /
  (panel_membership(emb[["Stim48hr"]], "Stim48hr") | panel_hubs(emb[["Stim48hr"]], "Stim48hr")) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 13, face = "bold")))

outfile <- file.path(script_dir, "..", "..", "plots", "knn_immune_target_manuscript.png")
# ragg's AGG device renders the Unicode glyphs the base PNG device drops;
# fall back to cairo if ragg is unavailable.
dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "cairo"
ggsave(outfile, fig, width = 9.8, height = 9.0, dpi = 300, bg = "white", device = dev)
cat("wrote", outfile, "\n")
