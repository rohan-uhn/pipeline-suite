#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
})

# CONFIG

CNVKIT_DIR <- commandArgs(trailingOnly = TRUE)[1]
OUT_TSV    <- commandArgs(trailingOnly = TRUE)[2]

if (is.na(CNVKIT_DIR) || is.na(OUT_TSV)) {
  stop("Usage: cnvkit_tp53_weighted_log2.R <cnvkit_dir> <output.tsv>")
}

TP53_CHR   <- "chr17"
TP53_START <- 7668421L
TP53_END   <- 7687490L

# HELPERS

weighted_mean <- function(x, w) {
  if (length(x) == 0 || sum(w, na.rm = TRUE) == 0) return(NA_real_)
  sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

# DISCOVER SAMPLES

sample_dirs <- list.dirs(CNVKIT_DIR, recursive = FALSE, full.names = TRUE)

# Keep only dirs that contain a .cnr file
sample_dirs <- sample_dirs[
  file.exists(file.path(sample_dirs,
                         paste0(basename(sample_dirs), "_realigned_recalibrated.cnr")))
]

if (length(sample_dirs) == 0) {
  stop("No CNVkit sample directories with .cnr files found in ", CNVKIT_DIR)
}

message("Found ", length(sample_dirs), " CNVkit samples")

# PROCESS EACH SAMPLE

rows <- vector("list", length(sample_dirs))

for (i in seq_along(sample_dirs)) {

  samp_dir <- sample_dirs[i]
  sample   <- basename(samp_dir)

  cnr_file <- file.path(
    samp_dir,
    paste0(sample, "_realigned_recalibrated.cnr")
  )

  if (!file.exists(cnr_file)) next

  dt <- fread(cnr_file)

  # Ensure required columns exist
  req <- c("chromosome", "start", "end", "log2", "weight")
  if (!all(req %in% colnames(dt))) {
    warning("Skipping ", sample, ": missing required CNVkit columns")
    next
  }

  # Overlap TP53
  tp53 <- dt[
    chromosome == TP53_CHR &
    start < TP53_END &
    end   > TP53_START
  ]

  if (nrow(tp53) == 0) {
    rows[[i]] <- data.table(
      sample           = sample,
      tp53_log2_weighted = NA_real_,
      n_bins           = 0L,
      weight_sum       = NA_real_,
      tp53_log2_median = NA_real_,
      tp53_log2_IQR    = NA_real_
    )
    next
  }

  rows[[i]] <- data.table(
    sample             = sample,
    tp53_log2_weighted = weighted_mean(tp53$log2, tp53$weight),
    n_bins             = nrow(tp53),
    weight_sum         = sum(tp53$weight, na.rm = TRUE),
    tp53_log2_median   = median(tp53$log2, na.rm = TRUE),
    tp53_log2_IQR      = IQR(tp53$log2, na.rm = TRUE)
  )
}

# OUTPUT

out <- rbindlist(rows, fill = TRUE)
setorder(out, tp53_log2_weighted)

fwrite(out, OUT_TSV, sep = "\t")

message("Wrote TP53 CNV summary to:")
message("  ", OUT_TSV)
