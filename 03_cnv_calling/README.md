# 03 – CNV Calling

This module performs copy number variant (CNV) analysis using a combination of:
1. **CNVkit (tumor-only, target-only mode)**
2. **Custom cohort-based CNV calling**
3. **Consensus CNV calling**

This directory contains **standalone scripts** and does **not** rely on the
PughLab pipeline framework beyond inputs generated in **02 – Preprocessing**.

> **Note:**  
> The current implementation focuses on **TP53**  
> (`chr17:7668421–7687490`, hg38), but some changes can generalize
> this workflow to other genes or regions.

---

## Inputs

This module requires outputs from **02 – Preprocessing**:

- Final realigned and recalibrated BAM files
- `gatk_bam_config_<timestamp>.yaml`
- Target BED files produced during preparation

---

## Software Environment

CNVkit analysis must be run inside the **CNVkit conda environment**:

```bash
conda activate cnvkit_env
````

Environment definition:

```
environments/cnvkit_env.yml
```

---

## Overview of CNVkit Workflow

The CNVkit portion of this module consists of four steps:

1. Create a **clean target BED** and **flat CNVkit reference**
2. Run CNVkit **per sample** in tumor-only, target-only mode
3. Verify CNVkit outputs
4. Aggregate CNVkit results across samples for the **TP53 region**

---

## 3.1 Build Clean Targets and Flat CNVkit Reference

This step:
- Cleans and validates the target BED file
- Generates CNVkit target bins
- Builds a **flat reference** (no matched normals)

### Script

```bash
bash make_cnvkit_reference_cleaned.sh \
  <gatk_bam_config.yaml> \
  <raw_targets.bed> \
  <output_directory> \
  <reference_genome.fa>
```

The script performs:

- BED cleaning (removal of malformed and duplicate intervals)
- BED sorting using reference `.fai`
- Target bin creation via `cnvkit.py target`
- Flat reference creation via `cnvkit.py reference`

### Outputs

```text
clean_targets.bed
targets.target.bed
flat_reference.cnn
```

These files are required for all downstream CNVkit runs.

---

## 3.2 Run CNVkit Per Sample (Tumor-Only)

CNVkit is run **independently per sample** using the BAM paths defined in
the GATK BAM configuration YAML.

### Script

```bash
bash run_cnvkit_per_sample.sh \
  <gatk_bam_config.yaml> \
  <targets.target.bed> \
  <flat_reference.cnn> \
  <output_directory>
```

### Key Characteristics

- Tumor-only mode
- Target-only analysis
- CNVkit scatter and diagram outputs enabled
- One SLURM job submitted per sample

### SLURM Configuration (Default)

```text
Time: 12 hours
CPUs: 4
Memory: 16 GB
```

Each sample is written to its own output directory:

```text
cnvkit/
├── SAMPLE1/
│   ├── SAMPLE1_realigned_recalibrated.cnr
│   ├── *.cns
│   ├── scatter.pdf
│   └── logs/
├── SAMPLE2/
└── ...
```

The `.cnr` files are the primary inputs for downstream CNV summarization.

---

## 3.3 Validate CNVkit Outputs

After CNVkit jobs complete, verify that all samples produced `.cnr` files.

### Script

```bash
bash check_targetcoverage.sh <cnvkit_output_directory>
```

### Output

The script reports:
- Samples with valid CNVkit outputs
- Samples missing `.cnr` files

This step should be completed **before proceeding** to cohort-level analysis.

---

## 3.4 Aggregate CNVkit Results for TP53

This step extracts CNVkit log2 ratios overlapping the **TP53 locus**
and summarizes them per sample using weighted statistics.

### TP53 Coordinates (hg38)

```text
chr17:7668421–7687490
```

---

### Script

```bash
Rscript cnvkit_tp53_weighted_log2.R \
  <cnvkit_output_directory> \
  <output.tsv>
```

---

### Metrics Computed Per Sample
- Weighted mean log2 ratio
- Median log2 ratio
- Interquartile range (IQR)
- Number of CNV bins overlapping TP53
- Total CNVkit weight across bins

---

### Output

```text
sample
tp53_log2_weighted
n_bins
weight_sum
tp53_log2_median
tp53_log2_IQR
```

This file forms the **CNVkit component** of the downstream
consensus CNV calling workflow.

---
