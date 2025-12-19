# 02 – Preprocessing

This module performs read-level quality control, alignment, and post-alignment
processing. It produces **final, analysis-ready BAM files** and associated
quality control metrics that are required for downstream CNV and SNV calling.

---

## Overview

The preprocessing stage consists of two major steps:

1. **Preprocessing pipeline**
   - FASTQC
   - Alignment with BWA-MEM
   - GATK realignment and base quality score recalibration (BQSR)

2. **BAM-level quality control**
   - Coverage and sequencing metrics using BAMQC

These steps are executed using selected components of the PughLab
DNA-seq pipeline framework.

---

## Inputs

This module requires the following inputs generated in **01 – Preparation**:

- FASTQ YAML configuration file  
  (`output_fastq_config.yaml`)
- Pipeline configuration YAML  
  (`configs/dna_pipeline_config.yaml`)
- Formatted target BED and interval files

---

## Output Directory Setup

Create an output directory for the cohort before running preprocessing:

```bash
mkdir <output_directory>
````

This directory will contain all preprocessing outputs, including FASTQC,
alignment results, and BAMQC metrics.

---

## 1. Preprocessing Pipeline (FASTQC, Alignment, GATK)

This step performs:

- Read-level QC using FASTQC
- Alignment using BWA-MEM
- GATK realignment and base quality score recalibration

---

### 1.1 Run Preprocessing

```bash
module load perl

perl <path_to_pipeline_suite>/pughlab_dnaseq_pipeline.pl \
  -t <path_to_pipeline_suite>/configs/dna_pipeline_config.yaml \
  -d <path_to_output_fastq_config.yaml> \
  --preprocessing \
  -c slurm \
  --remove
```

#### Notes

- Jobs are submitted to the HPC scheduler (SLURM)
- The pipeline can be safely re-run if jobs fail
- Execution will resume from the last completed step

---

### 1.2 Monitoring Job Status

After submission, inspect SLURM log files in the `log` folders created:

```text
slurm_job_metrics_*.out
```

If jobs fail due to:

- Out-of-memory errors
- Node failures
- Time limits
- Incorrect inputs

identify the affected samples/problem and re-run the pipeline. The workflow
will resume from the last successful step.

---

## 2. FASTQC Outputs

FASTQC results are written to a `fastqc/` directory within the output folder.

### Structure

```text
fastqc/
├── SAMPLE1/
├── SAMPLE2/
└── ...
```

- Summary reports provide an overview of read quality across the cohort
- Sample-level directories contain detailed FASTQC outputs

These results are useful for identifying:

- Low-quality samples
- Adapter contamination
- Sequencing anomalies

---

## 3. Alignment and GATK Outputs

### Final BAM Files

After preprocessing, **final analysis-ready BAM files** are located in:

```text
GATK/
├── SAMPLE1/
│   └── SAMPLE1_realigned_recalibrated.bam
├── SAMPLE2/
│   └── SAMPLE2_realigned_recalibrated.bam
├── gatk_bam_config_<DATE>_<TIME>.yaml
└── ...
```

These BAM files are the primary inputs for **CNV and SNV calling**.

---

### BAM Configuration YAML

A BAM configuration YAML file is generated automatically, with a timestamped name:

```text
gatk_bam_config_<DATE>_<TIME>.yaml
```

#### Example Structure

```yaml
---
SAMPLE1:
  tumour:
    SAMPLE1: /path/to/GATK/SAMPLE1/SAMPLE1_realigned_recalibrated.bam
SAMPLE2:
  tumour:
    SAMPLE2: /path/to/GATK/SAMPLE2/SAMPLE2_realigned_recalibrated.bam
```

This file maps each sample to its final BAM file location and is a **critical input**
for downstream CNV and SNV workflows.

---

## 4. BAM-Level Quality Control (BAMQC)

BAMQC is run separately after preprocessing to generate coverage and
sequencing metrics across the cohort.

---

### 4.1 Run BAMQC

```bash
module load perl

perl <path_to_pipeline_suite>/pughlab_dnaseq_pipeline.pl \
  -t <path_to_dna_pipeline_config.yaml> \
  -d <path_to_gatk_bam_config.yaml> \
  --qc \
  -c slurm \
  --remove
```

---

## 5. BAMQC Outputs

BAMQC results are written to a `BAMQC/` directory containing two main components:

```text
BAMQC/
├── Coverage/
├── SequenceMetrics/
```

---

### Coverage Metrics

The `Coverage/` directory contains:

* **Cohort-level summary file**, e.g.:

```text
<date>_<project>_Coverage_summary.tsv
```

This file reports:

- Mean coverage per sample
- Coverage distribution statistics across target regions
- **Sample-level coverage files**, including:
  - `*_DepthOfCoverage.read_group_interval_summary`
  - `*_DepthOfCoverage.sample_interval_summary`

These files are used in:

- Sample-level QC
- Coverage-based filtering
- CNV calling and interpretation (Part 03)

---

### Sequence Metrics

The `SequenceMetrics/` directory contains additional sequencing statistics,
including:

- Insert size metrics
- Read distribution summaries

---

## Outputs Summary

After completing this module, the following outputs are required for
downstream analysis:

- Final realigned and recalibrated BAM files
- `gatk_bam_config_<timestamp>.yaml`
- BAMQC coverage summary files
- Sample-level depth of coverage metrics

---

## Next Step

Proceed to:

→ **`03_cnv_calling/`**
for cohort-based CNV detection and summarization.

