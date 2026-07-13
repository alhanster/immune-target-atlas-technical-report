#!/usr/bin/env Rscript
# make_trial_fig.R
#
# Two-panel clinical-development validation figure (ggplot2 / patchwork)
# reproducing the figure that make_trial_fig.py draws. This is a pure RENDERER:
# it reads the two already-written validation data files and does NOT re-pull
# from Open Targets (run make_trial_fig.py to regenerate those files). Theme
# mirrors Part 2's make_regulator_burden_figure.R.
#
#   Panel a  top-25 drug programs -- for each of the top 25 unlabeled PU
#                                    nominations, the highest clinical stage of
#                                    any drug against that target, with a red
#                                    star on targets pursued for an immune
#                                    indication.
#   Panel b  enrichment vs control -- fraction of genes with a drug in trials,
#                                    top-25 nominations vs random lower-ranked
#                                    controls, for any and immune indications,
#                                    with Fisher exact p-values.
#
# Inputs (path resolved relative to this script; written by make_trial_fig.py):
#   trial_validation_top25.csv     pu_rank, gene, max_clinical_stage,
#                                  immune_max_stage, ...
#   trial_validation_counts.json   top_n, top_any_drug, top_immune, ctrl_n,
#                                  ctrl_any_drug, ctrl_immune
#
# Output (written next to this script):
#   trial_manuscript.png   (300 dpi; does not overwrite the Python-made
#                           trial_fig.png)
#
# Usage:  Rscript make_trial_fig.R

# Force a UTF-8 locale so Unicode glyphs in labels (× − ★) render instead of
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
  p <- normalizePath(file.path(script_dir, "..", "..", "data", "derived", f), mustWork = FALSE)
  if (!file.exists(p))
    stop("Missing ", p,
         "\nRun `python make_trial_fig.py` (or `--fetch`) first to write the ",
         "validation data files.")
  p
}

# --- Shared style (mirrors Part 2's make_regulator_burden_figure.R) ---------
col_top   <- "#2a7ab0"   # top nominations (blue)
col_ctrl  <- "#b0b0b0"   # random controls (grey)
col_star  <- "#d1495b"   # immune-indication star (red)

# ordinal phase -> numeric x position on the panel-a axis (from the .py).
# The "" (no-drug) case is handled by the NA fall-through below, not a named
# entry -- R forbids an empty ("") name in c().
STAGE_NUM <- c("Approved" = 4, "Phase IV" = 4, "Phase III" = 3, "Phase II/III" = 2.5,
               "Phase II" = 2, "Phase I/II" = 1.5, "Phase I" = 1, "Early Phase I" = 0.5)
# single-hue ordinal ramp: deeper blue = later stage (from the .py)
STAGE_COL <- c("Approved" = "#08519c", "Phase IV" = "#08519c", "Phase III" = "#2a7ab0",
               "Phase II/III" = "#4a97c9", "Phase II" = "#74add1", "Phase I/II" = "#a6bddb",
               "Phase I" = "#c6dbef")

# Font with full Unicode coverage so × − ★ render (macOS "sans" drops them).
font_family <- "Arial Unicode MS"
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

base_theme <- theme_classic(base_size = 9) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = 9),
        axis.text   = element_text(size = 7),
        plot.title  = element_text(size = 9, hjust = 0),
        plot.margin = margin(6, 10, 6, 10))

# ===================== Panel a: top-25 highest clinical stage =====================
top <- read.csv(need("trial_validation_top25.csv"), stringsAsFactors = FALSE)
top[is.na(top)] <- ""                                     # mirror .fillna("")
top <- top[order(top$pu_rank), ]                          # best (lowest rank) first

stage    <- as.character(top$max_clinical_stage)
max_num  <- unname(STAGE_NUM[stage]); max_num[is.na(max_num)] <- -1
ims      <- trimws(as.character(top$immune_max_stage))
top$y        <- rev(seq_len(nrow(top)))                   # top nomination at the top
top$plot_val <- ifelse(max_num < 0, 0.10, max_num)        # stub bar for "no known drug"
top$fill     <- ifelse(max_num < 0, "#f0f0f0",
                       ifelse(is.na(STAGE_COL[stage]), "#e3e3e3", STAGE_COL[stage]))
top$is_immune<- nzchar(ims) & ims != "nan"
top$star_x   <- pmax(max_num, 0) + 0.16

panelA <- ggplot(top, aes(x = plot_val, y = y)) +
  geom_col(aes(fill = fill), orientation = "y", width = 0.68,
           colour = "white", linewidth = 0.3) +
  geom_text(data = top[top$is_immune, ], aes(x = star_x, y = y),
            label = "★", colour = col_star, size = 2.9) +           # immune star
  # in-panel legend in the empty lower-mid whitespace (mirrors the .py)
  annotate("text", x = 2.55, y = 4.0, label = "★",
           colour = col_star, size = 2.9) +
  annotate("text", x = 2.72, y = 4.0, hjust = 0, size = 2.4,
           label = "= in trials for an\n   immune indication") +
  scale_fill_identity() +
  scale_x_continuous(breaks = 0:4, labels = c("none", "Ph I", "Ph II", "Ph III", "Appr"),
                     limits = c(0, 4.7), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = top$y, labels = top$gene,
                     expand = expansion(add = 0.6)) +
  labs(x = "Highest clinical stage of a drug against the target", y = NULL,
       title = "Top 25 novel nominations: existing drug programs", tag = "a") +
  base_theme +
  theme(axis.text.y = element_text(size = 6.6, face = "italic"),
        plot.tag    = element_text(size = 13))

# ==================== Panel b: enrichment vs control ====================
counts <- jsonlite::fromJSON(need("trial_validation_counts.json"))
tn <- counts$top_n;  cn <- counts$ctrl_n
t_any <- counts$top_any_drug; c_any <- counts$ctrl_any_drug
t_imm <- counts$top_immune;   c_imm <- counts$ctrl_immune
p_any <- fisher.test(matrix(c(t_any, tn - t_any, c_any, cn - c_any), nrow = 2, byrow = TRUE))$p.value
p_imm <- fisher.test(matrix(c(t_imm, tn - t_imm, c_imm, cn - c_imm), nrow = 2, byrow = TRUE))$p.value

cats     <- c("Any drug\nin trials", "Immune-indication\ntrial drug")
grp_top  <- sprintf("Top %d nominations", tn)
grp_ctrl <- sprintf("Random controls (n=%d)", cn)
# numeric x for the two categories; dodge the two groups by a fixed offset so the
# numeric p-value / label annotations share one continuous x scale.
barsB <- rbind(
  data.frame(x = c(1, 2) - 0.18, pct = c(100 * t_any / tn, 100 * t_imm / tn),
             group = grp_top,  lab_col = "black"),
  data.frame(x = c(1, 2) + 0.18, pct = c(100 * c_any / cn, 100 * c_imm / cn),
             group = grp_ctrl, lab_col = "#666666")
)
barsB$group <- factor(barsB$group, levels = c(grp_top, grp_ctrl))

panelB <- ggplot(barsB, aes(x = x, y = pct, fill = group)) +
  geom_col(width = 0.36) +
  geom_text(aes(label = sprintf("%.0f%%", pct), colour = lab_col),
            vjust = -0.5, size = 2.7, show.legend = FALSE) +
  annotate("text", x = 1, y = 26.5, label = sprintf("p=%.3f", p_any),
           size = 2.3, colour = "#08519c") +
  annotate("text", x = 2, y = 24.5, label = sprintf("p=%.3f", p_imm),
           size = 2.3, colour = "#08519c") +
  scale_fill_manual(values = setNames(c(col_top, col_ctrl), c(grp_top, grp_ctrl)),
                    name = NULL) +
  scale_colour_identity() +
  scale_x_continuous(breaks = c(1, 2), labels = cats, limits = c(0.5, 2.5)) +
  scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% of genes",
       title = "Nominations enriched for\nactive drug programs", tag = "b") +
  base_theme +
  theme(axis.text.x            = element_text(size = 7),
        legend.position        = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification   = c(1, 1),
        legend.text            = element_text(size = 6.2, family = font_family),
        legend.key.size        = unit(9, "pt"),
        plot.tag               = element_text(size = 13))

# ============================== Compose + save ==============================
fig <- panelA + panelB + plot_layout(widths = c(2.1, 1))

outfile <- file.path(script_dir, "..", "..", "plots", "trial_manuscript.png")
# ragg's AGG device renders the Unicode glyphs the base PNG device drops;
# fall back to cairo if ragg is unavailable.
dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "cairo"
ggsave(outfile, fig, width = 9.8, height = 5.2, dpi = 300, bg = "white", device = dev)
cat("wrote", outfile, "\n")
