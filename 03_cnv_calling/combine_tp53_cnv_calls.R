#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: combine_tp53_cnv_calls.R <cnvkit_tp53.tsv> <custom_tp53.tsv> <output.tsv>")
}

cnvkit_file <- args[1]
custom_file <- args[2]
out_file    <- args[3]

# Read inputs
cnvkit <- fread(cnvkit_file)
custom <- fread(custom_file)

# Standardize / select columns
cnvkit <- cnvkit[, .(
  sample,
  cnvkit_tp53_log2   = tp53_log2_weighted,
  cnvkit_n_bins      = n_bins,
  cnvkit_weight_sum = weight_sum,
  cnvkit_log2_median = tp53_log2_median,
  cnvkit_log2_IQR    = tp53_log2_IQR
)]

custom <- custom[, .(
  sample,
  custom_tp53_log2     = tp53_gene_log2,
  custom_z_gene        = z_gene,
  custom_log2_IQR      = tp53_log2_IQR,
  custom_n_intervals   = n_tp53_intervals
)]

# Recompute calls
cnvkit[, cnvkit_call :=
  ifelse(cnvkit_tp53_log2 < -0.30, "Loss", "Neutral")
]

custom[, custom_call :=
  ifelse(custom_tp53_log2 < -0.18, "Loss", "Neutral")
]

# Merge
merged <- merge(
  cnvkit,
  custom,
  by = "sample",
  all = TRUE
)

# Missingness annotation
merged[, missing_method :=
  ifelse(
    is.na(cnvkit_tp53_log2) & is.na(custom_tp53_log2), "both_missing",
    ifelse(
      is.na(cnvkit_tp53_log2), "cnvkit_missing",
      ifelse(is.na(custom_tp53_log2), "custom_missing", "none")
    )
  )
]

# Agreement
merged[, call_agreement :=
  ifelse(
    missing_method != "none",
    NA_character_,
    ifelse(cnvkit_call == custom_call, "MATCH", "DISCORDANT")
  )
]

# Write output
fwrite(merged, out_file, sep = "\t", quote = FALSE, na = "NA")

# Summary stats
cat("\n=== TP53 CNV comparison summary ===\n")
cat("Total samples:", nrow(merged), "\n")
cat("Both present:", sum(merged$missing_method == "none"), "\n")
cat("CNVkit missing:", sum(merged$missing_method == "cnvkit_missing"), "\n")
cat("Custom missing:", sum(merged$missing_method == "custom_missing"), "\n")
cat("Both missing:", sum(merged$missing_method == "both_missing"), "\n\n")

cat("Call agreement (where both present):\n")
print(table(merged$call_agreement, useNA = "ifany"))

cat("\nDone.\n")
