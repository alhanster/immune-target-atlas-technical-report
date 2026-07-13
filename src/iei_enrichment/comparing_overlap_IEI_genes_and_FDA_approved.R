#!/usr/bin/env Rscript
# comparing_overlap_IEI_genes_and_FDA_approved.R
#
# Compares the Inborn Errors of Immunity (IEI) gene list against the targets of
# FDA-approved, immune-indicated therapeutics, and reports the overlap at both
# the gene level and the drug (therapeutic) level.
#
# Inputs (paths relative to the repo root):
#   Data/Gene List/master_gene_list.csv
#       IEI-associated genes (column `Gene`, `Inheritance`). HGNC-normalized.
#   Data/Drug List/fda_approved_target_genes.txt
#       Unique target genes of FDA-approved drugs with an immune indication
#       (one HGNC symbol per line).
#   Data/Drug List/Other/immune_system_drugs_fda_approved.csv
#       Drug-level table (Open Targets derived, maxClinicalStage == APPROVAL).
#       Used only for the drug-level (therapeutic) count; the immune-indicated
#       subset is rows where `immune_indication == 1`.
#
# Outputs:
#   Printed summary table, plus:
#     gene_overlap_summary.csv   (the summary metrics)
#     overlapping_genes.csv      (genes present in both lists)
#
# Usage:  Rscript comparing_overlap_IEI_genes_and_FDA_approved.R

# --- Resolve paths relative to this script's location -----------------------
# This script lives in pre-lim-result/; its parent is the repo root. Resolving
# from the script location lets it run from any working directory.
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_root   <- dirname(script_dir)

# --- Input file paths -------------------------------------------------------
iei_path   <- file.path(repo_root, "Data", "Gene List", "master_gene_list.csv")
fda_path   <- file.path(repo_root, "Data", "Drug List", "fda_approved_target_genes.txt")
drug_path  <- file.path(repo_root, "Data", "Drug List", "Other",
                        "immune_system_drugs_fda_approved.csv")

# --- Load gene sets ---------------------------------------------------------
# The master list has been HGNC-normalized (no BOM / stray whitespace), so a
# plain read.csv is sufficient. We still strip any residual BOM defensively.
iei_lines <- readLines(iei_path, warn = FALSE)
iei_lines[1] <- sub("^﻿", "", iei_lines[1])
iei_tbl <- read.csv(text = iei_lines, stringsAsFactors = FALSE, check.names = FALSE)
iei_genes <- unique(trimws(iei_tbl$Gene))
iei_genes <- iei_genes[iei_genes != ""]

fda_genes <- unique(trimws(readLines(fda_path, warn = FALSE)))
fda_genes <- fda_genes[fda_genes != ""]

# --- Gene-level overlap -----------------------------------------------------
in_both     <- sort(intersect(iei_genes, fda_genes))
iei_only    <- sort(setdiff(iei_genes, fda_genes))
fda_only    <- sort(setdiff(fda_genes, iei_genes))
union_genes <- union(iei_genes, fda_genes)

pct <- function(n, d) if (d == 0) NA_real_ else round(100 * n / d, 1)

# --- Drug-level (therapeutic) overlap ---------------------------------------
# Of the FDA-approved, immune-indicated drugs, how many target >= 1 IEI gene?
# A single target gene is hit by many drugs, so this is a genuinely different
# denominator from the gene-level figures above.
drug_pct <- NA_real_; n_drugs <- NA; n_drugs_iei <- NA
if (file.exists(drug_path)) {
  d <- read.csv(drug_path, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[d$immune_indication == 1, ]
  drug <- d[["Drug Name"]]; gene <- d[["Target (Gene)"]]
  n_drugs     <- length(unique(drug))
  n_drugs_iei <- length(unique(drug[gene %in% iei_genes]))
  drug_pct    <- pct(n_drugs_iei, n_drugs)
}

# --- Summary table ----------------------------------------------------------
summary_tbl <- data.frame(
  Metric = c(
    "IEI genes (unique)",
    "FDA immune-drug target genes (unique)",
    "Gene overlap (in both)",
    "IEI only",
    "FDA only",
    "Union (either list)",
    "% of IEI genes that are FDA immune-drug targets",
    "% of FDA target genes that are IEI genes",
    "Jaccard similarity (%)",
    "FDA-approved immune-indicated drugs (unique)",
    "Drugs targeting >=1 IEI gene",
    "% of therapeutics targeting >=1 IEI gene"
  ),
  Value = c(
    length(iei_genes),
    length(fda_genes),
    length(in_both),
    length(iei_only),
    length(fda_only),
    length(union_genes),
    pct(length(in_both), length(iei_genes)),   # 14.3%
    pct(length(in_both), length(fda_genes)),   # 10.0%
    pct(length(in_both), length(union_genes)),
    n_drugs,
    n_drugs_iei,
    drug_pct                                    # 17.9%
  ),
  stringsAsFactors = FALSE
)

# --- Report -----------------------------------------------------------------
cat("Overlap: IEI genes vs. FDA-approved immune-drug targets\n")
cat("=======================================================\n")
print(summary_tbl, row.names = FALSE, right = FALSE)

cat("\nNote on units: the gene-level figures (10.0% / 14.3%) count GENES;\n")
cat("the 17.9% figure counts DRUGS (therapeutics). They are not interchangeable.\n")

cat("\nFirst overlapping genes:\n")
cat(paste(head(in_both, 20), collapse = ", "),
    if (length(in_both) > 20) ", ..." else "", "\n", sep = "")

# --- Write outputs (next to this script) ------------------------------------
write.csv(summary_tbl, file.path(script_dir, "gene_overlap_summary.csv"), row.names = FALSE)
write.csv(data.frame(Gene = in_both),
          file.path(script_dir, "overlapping_genes.csv"), row.names = FALSE)
cat("\nWrote gene_overlap_summary.csv and overlapping_genes.csv to", script_dir, "\n")
