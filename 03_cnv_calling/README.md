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

## 3.1 Overview of CNVkit Workflow

The CNVkit portion of this module consists of four steps:

1. Create a **clean target BED** and **flat CNVkit reference**
2. Run CNVkit **per sample** in tumor-only, target-only mode
3. Verify CNVkit outputs
4. Aggregate CNVkit results across samples for the **TP53 region**

---

## 3.1.1 Build Clean Targets and Flat CNVkit Reference

This step:
- Cleans and validates the target BED file
- Generates CNVkit target bins
- Builds a **flat reference** (no matched normals)

### Script

```bash
bash <path_to_make_cnvkit_reference_cleaned.sh> \
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

## 3.1.2 Run CNVkit Per Sample (Tumor-Only)

CNVkit is run **independently per sample** using the BAM paths defined in
the GATK BAM configuration YAML.

### Script

```bash
bash <path_to_run_cnvkit_per_sample.sh> \
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

## 3.1.3 Validate CNVkit Outputs

After CNVkit jobs complete, verify that all samples produced `.cnr` files.

### Script

```bash
bash <path_to_check_targetcoverage.sh> <cnvkit_output_directory>
```

### Output

The script reports:
- Samples with valid CNVkit outputs
- Samples missing `.cnr` files

This step should be completed **before proceeding** to cohort-level analysis.

---

## 3.1.4 Aggregate CNVkit Results for TP53

This step extracts CNVkit log2 ratios overlapping the **TP53 locus**
and summarizes them per sample using weighted statistics.

### TP53 Coordinates (hg38)

```text
chr17:7668421–7687490
```

---

### Script

```bash
module load R/3.6.1

Rscript <path_to_cnvkit_tp53_weighted_log2.R> \
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

## 3.2 Custom Cohort-Based CNV Calling (TP53)

This step performs **custom cohort-based CNV calling** for a specified genomic
region using **BAMQC coverage outputs**, independent of CNVkit.

The goal of this workflow is to derive a **gene-level copy number signal**
by leveraging:
- Within-sample coverage normalization
- Cohort-level reference construction
- Robust interval aggregation

This approach provides an **independent CNV signal** that complements CNVkit
and is later integrated during **consensus CNV calling**.

> **Current scope:**  
> This implementation is configured for **TP53**  
> (`chr17:7668421–7687490`, hg38), but can be generalized to other regions
> with minimal modification.

---

## Rationale

CNVkit operates on read-depth segmentation and bin-level modeling.
In contrast, this custom workflow:
- Uses **target-level coverage directly from BAMQC**
- Normalizes coverage **within each sample**
- Builds a **reference from the cohort itself**
- Is robust to tumor-only sequencing data
- Enables explicit **sample-level QC and exclusion**

This provides an orthogonal view of copy number state that is less sensitive
to segmentation artifacts and more interpretable at the gene level.

---

## Inputs

This step requires:

- BAMQC outputs from **02 – Preprocessing**
- Specifically, per-sample files:
```
<SAMPLE>_DepthOfCoverage.sample_interval_summary
```
- A user-defined list of samples to include in the cohort

---

## 3.2.1 Sample Selection and Filtering

### Purpose

Samples with very low coverage can bias cohort-based CNV estimates.
To avoid this, samples are **explicitly filtered prior to analysis**.

### Procedure

1. Inspect the BAMQC coverage summary generated in Part 2:
```
BAMQC/Coverage/<date>_<project>_Coverage_summary.tsv
```

2. Identify samples meeting a desired coverage threshold
(e.g. mean coverage > 50×)
3. Create a plain text file containing **one sample ID per line**:

```text
SAMPLE1
SAMPLE2
SAMPLE3
````

These sample IDs must correspond to BAMQC subdirectories containing:

```text
<SAMPLE>_DepthOfCoverage.sample_interval_summary
```

This file controls which samples are included in the cohort analysis.

---

## 3.2.2 Wrapper Script Execution

A shell wrapper script is used to execute the cohort-based CNV workflow.

### Wrapper Usage

```bash
bash <path_to_run_tp53_cnv_wrapper.sh> \
  <sample_list.txt> \
  <bamqc_dir> \
  <outdir>
```

### Wrapper Behavior

- Loads an R environment
- Executes the cohort-based CNV R script
- Writes all outputs to the specified output directory
- Fails early if required inputs are missing

---

## 3.2.3 Cohort-Based CNV Algorithm

The R script (`tp53_cnv_from_bamqc_hg38_from_list.R`) performs the following steps.

---

### Step 1: Load BAMQC Interval Coverage

For each sample:
- Read per-target coverage from:
  ```
  <sample>_DepthOfCoverage.sample_interval_summary
  ```
- Extract:
  - Chromosome
  - Interval start and end
  - Mean coverage
- Compute interval length

These values form the raw coverage signal for CNV inference.

---

### Step 2: Within-Sample Normalization

To account for differences in sequencing depth and library size, coverage
values are normalized **within each sample**.

For a given sample *i*:
- A scaling factor is computed as the **median non-zero coverage**
- Each interval’s coverage is divided by this scaling factor

Conceptually:
```text
normalized_coverage = interval_coverage / median_sample_coverage
```

This ensures that downstream copy number estimates reflect **relative depth
changes**, rather than absolute sequencing depth.

---

### Step 3: Cohort Reference Construction

A cohort reference is constructed using normalized coverage values:

- Samples ending in `_T2` are excluded from reference construction
- For each interval across the cohort:
  - Compute the cohort median normalized coverage
  - Compute the median absolute deviation (MAD)
  - Exclude outlier intervals where deviation > 3 × MAD

The remaining intervals define a **robust cohort reference** representing
expected copy-neutral coverage.

---

### Step 4: Copy Ratio Calculation

For each sample and interval, a copy ratio (CR) is computed as:
```text
CR = normalized_coverage / reference_coverage
```

This is converted to log2 scale:
```text
log2CR = log2(CR)
```

A negative log2CR indicates relative copy loss, while positive values
indicate relative copy gain.

---

### Step 5: TP53 Interval Selection

Only intervals overlapping the TP53 locus:
```text
chr17:7668421–7687490
```

are retained for gene-level summarization.

---

### Step 6: Gene-Level Aggregation

TP53 interval-level values are collapsed to a single gene-level estimate
per sample.

This is done using a **length-weighted median** of interval-level log2 copy ratios:

- Longer intervals contribute more weight
- Outlier intervals have reduced influence
- The result is robust to noisy targets

Additional summary statistics are also reported:

- Median log2 copy ratio
- Interquartile range (IQR)
- Number of contributing intervals

---

### Step 7: Z-Score Normalization (Interpretation Only)

To contextualize gene-level copy ratios within the cohort, a **z-score**
is computed across samples.

For each sample *i*:

- Let (x_i) be the TP53 gene-level log2 copy ratio
- Let (μ) be the cohort median of ( x )
- Let (σ) be the cohort standard deviation of ( x )

The z-score is computed as:

```text
z_i = (x_i − μ) / σ
```

A provisional loss classification is assigned as:

```text
z ≤ −0.56 → Loss
z > −0.56 → Neutral
```

---

#### Important Note on Usage

> **Z-score–based loss calls are NOT currently used as the primary CNV call
> in downstream analyses.**

Z-scores are retained as:
- A secondary signal
- A diagnostic measure of extremeness
- A potential alternative calling strategy if cohort-wide distributions
  are considered more informative than absolute log2 thresholds

Final CNV calls are made during **consensus CNV calling** (Part 3.3).

---

## 3.2.4 Outputs

The pipeline produces the following outputs.

### Gene-Level Summary
```text
tp53_gene_copyratio.tsv
```

Contains gene-level log2 copy ratios, z-scores, and summary statistics.

---

### Per-Interval Copy Ratios
```text
per_interval_copyratio.tsv.gz
```

Contains normalized coverage and copy ratio values for all intervals.

---

### Coverage QC Metrics
```text
depth_qc.tsv
```

Provides per-sample coverage and variability summaries.

---

### Final Filtered Output
```text
final/tp53_gene_copyratio_final.tsv
```

This file is used directly in **consensus CNV calling** (Part 3.3).

## 3.3 Consensus CNV Calling

This step performs **consensus CNV calling** for TP53 by integrating results from:

1. **CNVkit (tumor-only, target-only)**
2. **Custom cohort-based CNV workflow (BAMQC-based)**

The goal is to produce a **robust, conservative TP53 loss call** that minimizes
false positives in tumor-only data.

> **Scope note:**  
> This step is designed to call **copy number loss only**.  
> Copy number gains are not currently assessed, but could be incorporated
> with additional thresholds and logic.

---

## Rationale

Tumor-only CNV calling is inherently noisy due to:
- Lack of matched normals
- Variable tumor purity
- Cohort-specific depth effects

To address this, consensus calling is used to:
- Combine **two independent CNV signals**
- Apply **conservatively calibrated thresholds**
- Reduce false-positive loss calls

This approach prioritizes **specificity over sensitivity**, particularly
important when interpreting tumor suppressor genes such as TP53.

---

## 3.3.1 Threshold Calibration Strategy

### Overview

The loss thresholds used in this step were **empirically calibrated** using
a separate benchmarking analysis involving:

- A **tumor+normal cohort**
- CNV calls generated using:
  - GATK CNV with a panel of normals (PoN)
  - CNVkit with a normals-based reference

The goal of calibration was to determine tumor-only thresholds such that:

- **No normal samples** were called as TP53 loss
- Tumor-only calls **closely matched** calls from normal-referenced workflows
- Agreement with GATK CNV was maximized

---

### Optimization Criteria

Thresholds were selected to:
- Eliminate false-positive loss calls in normals
- Maximize **F1 score** when compared against:
  - GATK CNV (PoN-based)
  - CNVkit (normal-reference-based)

As a result, these thresholds are **not arbitrary**, but tuned to mimic
gold-standard CNV calling as closely as possible using tumor-only data.

---

## 3.3.2 Applied Loss Thresholds

The following **log2 copy ratio thresholds** are applied:

| Method    | Metric                     | Loss Threshold |
|-----------|----------------------------|----------------|
| CNVkit   | Weighted TP53 log2 ratio    | `< -0.30`     |
| Custom   | Gene-level TP53 log2 ratio  | `< -0.18`     |

Only **loss vs neutral** states are considered.

---

## 3.3.3 Inputs

This step requires two TSV files generated previously:

- **CNVkit TP53 summary**
```
tp53_cnvkit_summary.tsv

```
- **Custom cohort-based TP53 summary**
```
tp53_gene_copyratio_final.tsv

````

Both files must contain a `sample` column.

---

## 3.3.4 Running Consensus CNV Calling

Consensus calling is performed using a single R script.

### Command

```bash
Rscript combine_tp53_cnv_calls.R \
<cnvkit_tp53.tsv> \
<custom_tp53.tsv> \
<output.tsv>
````

---

## 3.3.5 Consensus Calling Logic

The script performs the following steps.

---

### Step 1: Input Harmonization

Relevant fields are extracted and standardized from each input:

**From CNVkit:**
- Weighted TP53 log2 ratio
- Number of bins
- Total bin weight
- Median and IQR of log2 ratios

**From Custom Workflow:**
- TP53 gene-level log2 ratio
- Z-score (for reference)
- IQR
- Number of contributing intervals

---

### Step 2: Independent Loss Calling

Loss calls are recomputed internally using calibrated thresholds:

```text
CNVkit:  log2 < −0.30 → Loss
Custom:  log2 < −0.18 → Loss
```

This ensures reproducibility and avoids reliance on upstream annotations.

---

### Step 3: Merge and Missingness Annotation

Results from both methods are merged by sample ID.

Each sample is annotated with a **missingness category**:

- `none` – both methods present
- `cnvkit_missing`
- `custom_missing`
- `both_missing`

This allows transparent downstream filtering.

---

### Step 4: Call Agreement Assessment

For samples where **both methods are present**, agreement is assessed:

- `MATCH` – both methods agree
- `DISCORDANT` – methods disagree

Samples with missing data from either method are not forced into agreement.

---

## 3.3.6 Outputs

### Primary Output

```text
<output.tsv>
```

This file contains, per sample:

- CNVkit TP53 metrics and call
- Custom TP53 metrics and call
- Missingness status
- Call agreement annotation

This is the **final integrated CNV table** used for downstream analyses,
visualization, and reporting.

---

### Console Summary

Upon completion, the script prints:

- Total number of samples
- Number of samples with:
  - Both methods present
  - CNVkit missing
  - Custom missing
  - Both missing
- Agreement vs discordance counts

This provides a quick QC check of cohort-level consistency.

---

## Interpretation Notes

* **Loss calls are intentionally conservative**
* Z-scores from the custom workflow are retained for context but are **not
  enforced** at this stage

---

## Extensibility

With modification, this framework can be extended to:
- Copy number **gain** calling
- Additional genes or regions
- Alternative threshold optimization strategies
- Weighted consensus schemes

---

## Next Step

Proceed to:

→ **04 – SNV Calling**

This stage performs multi-caller SNV detection and filtering.
