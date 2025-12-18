# 01 – Preparation

This module performs all setup steps required before running the preprocessing
and variant calling workflows. It prepares target BED files, defines sample
metadata, and generates configuration files used by downstream pipeline stages.

---

## Overview

The preparation step consists of four main tasks:

1. Formatting and standardizing the target BED file  
2. Organizing raw FASTQ files and generating a sample information table  
3. Updating the main pipeline configuration file  
4. Generating a FASTQ YAML configuration file  

These steps must be completed **once per cohort** before running the pipeline.

---

## Assumptions & Directory Structure

### FASTQ File Requirements

- Input data must start from **paired-end FASTQ files**
- FASTQ files must follow Illumina-style naming:
  - `SAMPLE_ID_R1_001.fastq.gz`
  - `SAMPLE_ID_R2_001.fastq.gz`
- All FASTQ files for a cohort should be placed in a single directory

### Recommended Directory Layout

```
cohort/
├── raw/                      # FASTQ files
│   ├── SAMPLE1_R1_001.fastq.gz
│   ├── SAMPLE1_R2_001.fastq.gz
│   └── ...
├── targets.bed               # Original target BED
├── sample_info.txt
└── output_fastq_config.yaml
```

## 1. Target BED File Formatting
The pipeline requires a sorted, indexed BED file as well as
derived interval and padded BED files.

---

### 1.1 Sort the BED File
```bash
sort -k1,1 -k2,2n -k3,3n \
  <targets.bed> \
  > <targets.sorted.bed>
```

---

### 1.2 Format BED and Generate Interval Files
This step uses a helper script from the PuLab pipeline-suite to:

- Validate and format the BED file
- Generate a Picard-compatible interval list
- Create a padded BED file (±100 bp)
- Compress and index BED outputs

```bash
module load perl

perl <path_to_pipeline_suite>/scripts/format_intervals_bed.pl \
  -b <targets.sorted.bed> \
  -r <reference_genome_fasta>
```

#### Expected Outputs
```text
targets.sorted.bed
targets.sorted.bed.gz
targets.sorted.bed.gz.tbi
targets.sorted.interval_list
targets.sorted_padding100bp.bed
targets.sorted_padding100bp.bed.gz
targets.sorted_padding100bp.bed.gz.tbi
```

These files are required for preprocessing and CNV calling.

---

## 2. Sample Information File
A sample information file is required to define:
- Sample identifiers
- Sample type (tumour)

For tumor-only datasets, **all samples are labeled as `tumour`**.

---

### 2.1 Generate `sample_info.txt`

The following commands assume:

- All FASTQ files are located in the `raw/` directory
- FASTQ filenames follow the required naming convention

```bash
echo "Patient.ID Sample.ID Type" > header.txt

ls raw/*R1_001.fastq.gz \
  | sed 's#^raw/##' \
  | sed 's/_R1_001.fastq.gz$//' \
  | awk '{print $1, $1, "tumour"}' \
  >> header.txt

cat header.txt | sed 's/ \+/\t/g' > sample_info.txt
rm header.txt
```

---

### Example Output

```text
Patient.ID	Sample.ID	Type
SAMPLE1	SAMPLE1	tumour
SAMPLE2	SAMPLE2	tumour
SAMPLE3	SAMPLE3	tumour
```

---

## 3. Update Pipeline Configuration File

The main pipeline configuration YAML `configs/dna_pipeline_config.yaml` must be updated to reflect the cohort-specific
paths and project name.

At minimum, the following fields should be updated:

```yaml
project_name: <cohort_name>
output_dir: <path_to_output_directory>
targets_bed: <path_to_sorted_targets.bed>
```

> **Note:**
> This configuration file primarily affects the preprocessing stage.
> CNV and SNV workflows use standalone scripts and conda environments.

---

## 4. Generate FASTQ YAML Configuration

The FASTQ YAML file defines the mapping between sample IDs and FASTQ files
and is required for downstream pipeline execution.

---

### 4.1 Create FASTQ YAML

```bash
module load perl

perl <path_to_pipeline_suite>/scripts/create_fastq_yaml.pl \
  -d <path_to_raw_fastq_directory> \
  -o <output_fastq_config.yaml> \
  -t dna \
  -i <sample_info.txt>
```

### Example Output Structure

```yaml
---
SAMPLE1:
  SAMPLE1:
    libraries:
      SAMPLE1:
        runlanes:
          001:
            fastq:
              R1: /path/to/raw/SAMPLE1_R1_001.fastq.gz
              R2: /path/to/raw/SAMPLE1_R2_001.fastq.gz
```

This file is consumed by the preprocessing pipeline to locate FASTQ inputs.

---

## Outputs Summary

After completing this module, the following files should be available:
- Formatted and indexed target BED files
- `sample_info.txt`
- Updated pipeline configuration YAML
- FASTQ YAML configuration file

These outputs are required inputs for **02 – Preprocessing**.

---

## Next Step

Proceed to:

→ **`02_preprocessing/`**
for quality control, alignment, and base recalibration.
