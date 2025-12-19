#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(stringr)
})

# INPUTS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript tp53_cnv_from_bamqc_hg38_from_list.R sample_list.txt bamqc_dir outdir")
}

sample_file <- args[1]
bamqc_dir   <- args[2]
outdir      <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "final"), showWarnings = FALSE, recursive = TRUE)

# LOAD SAMPLE LIST
samples <- fread(sample_file, header = FALSE)[[1]]
samples <- samples[nzchar(samples)]
samples <- unique(samples)

if (length(samples) == 0) stop("Sample list is empty.")

message("Loaded ", length(samples), " samples.")

# CONSTANTS
TP53_CHR   <- "chr17"
TP53_START <- 7668421L
TP53_END   <- 7687490L
CUT_Z      <- -0.56   # z-score cutoff for loss

# READ BAMQC COVERAGE
rows <- list()

for (samp in samples) {

  f <- file.path(bamqc_dir, samp, paste0(samp, "_DepthOfCoverage.sample_interval_summary"))

  if (!file.exists(f)) {
    warning("Missing BAMQC file for sample: ", samp)
    next
  }

  dt <- suppressWarnings(fread(f, sep="\t", header=TRUE, fill=TRUE))

  if (!all(c("Target", "average_coverage") %in% names(dt))) {
    warning("Missing required columns in file for sample: ", samp)
    next
  }

  sub <- dt[, .(interval = Target, mean_cov = as.numeric(average_coverage))]
  sub <- sub[!is.na(mean_cov)]

  iv <- tstrsplit(sub$interval, "[:\\-]")
  sub[, chr   := iv[[1]] ]
  sub[, start := as.integer(iv[[2]])]
  sub[, end   := as.integer(iv[[3]])]
  sub[, len   := pmax(end - start, 1L)]
  sub[, sample := samp]

  rows[[samp]] <- sub[, .(sample, chr, start, end, len, mean_cov)]
}

x <- rbindlist(rows, use.names = TRUE, fill = TRUE)
if (nrow(x) == 0) stop("No interval coverage data read.")

# WITHIN-SAMPLE NORMALIZATION
sf <- x[, .(sf = median(mean_cov[mean_cov > 0], na.rm=TRUE)), by=sample]
x <- merge(x, sf, by="sample")
x[, norm := mean_cov / pmax(sf, 1e-9)]

# BUILD COHORT REFERENCE (EXCLUDE _T2)
base_samples <- samples[!grepl("_T2$", samples)]
message("Building reference using ", length(base_samples), " samples (excluding T2).")

ref0 <- x[sample %in% base_samples,
          .(med = median(norm, na.rm=TRUE),
            mad = mad(norm, na.rm=TRUE)),
          by = .(chr,start,end)]

x <- merge(x, ref0, by=c("chr","start","end"))

x[, keep := abs(norm - med) <= 3 * pmax(mad, 1e-6)]

ref <- x[keep == TRUE & sample %in% base_samples,
         .(ref_interval = median(norm, na.rm=TRUE)),
         by = .(chr,start,end)]

x <- merge(x[, !"med"], ref, by=c("chr","start","end"), all.x=TRUE)

# COPY RATIO COMPUTATION
x[, CR := norm / pmax(ref_interval, 1e-9)]
x[, log2CR := log2(pmax(CR, 1e-9))]

# TP53 INTERVALS
x[, is_tp53 := (chr == TP53_CHR & start < TP53_END & end > TP53_START)]

tp53_iv <- x[is_tp53 == TRUE]
if (nrow(tp53_iv) == 0) stop("No intervals overlap TP53.")

# TP53 GENE-LEVEL COLLAPSE
weighted.median <- function(z, w) {
  o <- order(z); z <- z[o]; w <- w[o] / sum(w)
  cw <- cumsum(w); z[which.max(cw >= 0.5)]
}

tp53_gene <- tp53_iv[, .(
  tp53_gene_log2   = weighted.median(log2CR, len),
  tp53_log2_median = median(log2CR),
  tp53_log2_IQR    = IQR(log2CR),
  n_tp53_intervals = .N
), by = sample]

# Z-SCORE NORMALIZATION
mu  <- median(tp53_gene$tp53_gene_log2, na.rm=TRUE)
sdv <- sd(tp53_gene$tp53_gene_log2, na.rm=TRUE)

tp53_gene[, z_gene :=
            if (is.finite(sdv) && sdv > 0)
              (tp53_gene_log2 - mu) / sdv
            else NA_real_]

tp53_gene[, z_call := ifelse(z_gene <= CUT_Z, "Loss", "Neutral")]

# OUTPUTS
fwrite(tp53_gene[order(z_gene)],
       file.path(outdir, "tp53_gene_copyratio.tsv"), sep="\t")

fwrite(x[order(sample, chr, start),
         .(sample, chr, start, end, len, mean_cov, norm, ref_interval,
           CR, log2CR, is_tp53)],
       file.path(outdir, "per_interval_copyratio.tsv.gz"), sep="\t")

qc <- x[, .(
  median_cov = median(mean_cov, na.rm=TRUE),
  log2CR_IQR = IQR(log2CR, na.rm=TRUE),
  n_targets = .N
), by = sample]

fwrite(qc, file.path(outdir, "depth_qc.tsv"), sep="\t")

final_tbl <- tp53_gene[, .(
  sample,
  tp53_gene_log2,
  z_gene,
  z_call,
  tp53_log2_IQR,
  n_tp53_intervals
)]

fwrite(final_tbl[order(z_gene)],
       file.path(outdir, "final", "tp53_gene_copyratio_final.tsv"),
       sep="\t")

message("DONE. Outputs written to:")
message(" - ", file.path(outdir, "tp53_gene_copyratio.tsv"))
message(" - ", file.path(outdir, "per_interval_copyratio.tsv.gz"))
message(" - ", file.path(outdir, "depth_qc.tsv"))
message(" - ", file.path(outdir, "final/tp53_gene_copyratio_final.tsv"))
