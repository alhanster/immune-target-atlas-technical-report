#!/usr/bin/env Rscript
# make_gwas_violin.R
#
# Violin plot comparing the GWAS credible-set gene score (immune system disorder)
# between approved drug targets and non-targets, over the full gene list.
# Styled after syngap1_olink/analysis (config.R + utils.R + VolcanoPlot_Compiled.R),
# self-contained: config and save_figure() inlined below.

# Force a UTF-8 locale so Unicode glyphs render instead of one dot per byte.
for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# =============================================================================
# Config (inlined; mirrors syngap1_olink config.R / utils.R)
# =============================================================================
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()

OUTPUT_DIR        <- file.path(script_dir, "..", "..", "plots")
DATA_DIR          <- file.path(script_dir, "..", "..", "data")
BASE_TEXT_SIZE    <- 7
FIGURE_DPI        <- 300
TARGET_COLOR      <- "#009E73"   # approved drug target (Okabe-Ito green)
NONTARGET_COLOR   <- "grey"      # non-target
colorblind_palette <- c(TARGET_COLOR, NONTARGET_COLOR)   # matched to factor order

# Save a figure to OUTPUT_DIR with consistent units and DPI.
save_figure <- function(plot, filename, width, height) {
  filepath <- file.path(OUTPUT_DIR, filename)
  ggsave(filepath, plot = plot, width = width, height = height, units = "mm", dpi = FIGURE_DPI)
  message("Saved: ", filepath)
}

# =============================================================================
# Load data
# =============================================================================
full_gene_list <- read.delim(file.path(DATA_DIR, "derived", "full_gene_list.tsv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
approved_targets <- readLines(file.path(DATA_DIR, "reference", "approved_target_genes.txt"))
approved_targets <- trimws(approved_targets)
approved_targets <- approved_targets[nzchar(approved_targets)]

plot_data <- full_gene_list %>%
  mutate(gwas_score = suppressWarnings(as.numeric(gwas_score))) %>%
  filter(!is.na(gwas_score)) %>%
  mutate(is_target = gene %in% approved_targets)

n_target    <- sum(plot_data$is_target)
n_nontarget <- sum(!plot_data$is_target)
lab_target    <- sprintf("approved\ndrug target\n(n = %s)", format(n_target, big.mark = ","))
lab_nontarget <- sprintf("Non-target\n(n = %s)",                format(n_nontarget, big.mark = ","))

plot_data <- plot_data %>%
  mutate(Target = factor(ifelse(is_target, lab_target, lab_nontarget),
                         levels = c(lab_target, lab_nontarget)))

# =============================================================================
# Violin plot: approved drug targets vs non-targets
# =============================================================================
plot_violin <- ggplot(plot_data, aes(x = Target, y = gwas_score, fill = Target)) +
  geom_violin(trim = TRUE, color = NA, alpha = 0.7) +
  labs(title = "", x = "Group", y = "GWAS score for MONDO:0005046 - immune system disorder") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = colorblind_palette) +
  theme_classic(base_size = BASE_TEXT_SIZE) +
  theme(legend.position = "none",
        legend.key.size = unit(3, "mm"),
        axis.title.x = element_blank())

# =============================================================================
# Save
# =============================================================================
save_figure(plot_violin, "figure_gwas_violin.png", width = 120, height = 90)
