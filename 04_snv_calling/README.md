# 04 – SNV Calling

This module performs **tumor-only single nucleotide variant (SNV) calling**
using five independent variant callers:

- MuTect
- MuTect2
- VarScan
- VarDict
- Strelka

Each caller is run independently on the same set of tumor-only BAM files.
The resulting calls are then:
1. Filtered
2. Annotated using **VEP** to generate MAF files
3. Combined across callers at the sample and cohort level
4. Used for downstream filtering, analysis, and visualization in R

> **Important:**  
> These workflows are executed as **standalone scripts** and do not use the
> PughLab pipeline framework, however they were inspired by their framework.

---

## Inputs

This module requires outputs from **02 – Preprocessing**:

- `gatk_bam_config_<timestamp>.yaml`
- Final realigned and recalibrated BAM files
- Target BED file generated in **01 – Preparation**
- Reference genome FASTA (`.fa`) and index (`.fai`)
- Panel of Normals (PoN), the [1000 Genomes hg38 PoN](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON) is recommeded

---

## Software Environments

Two conda environments are required:

- **`cnvkit_env`**
  - Used for some Python and post-processing steps
- **`vep_env`**
  - Required for VEP annotation and `vcf2maf.pl`

Install from the repository:

```bash
conda env create -f environments/cnvkit_env.yml
conda env create -f environments/vep_env.yml
````

Each script will activate the appropriate environment.

---

## General Workflow

For each caller, the workflow follows the same high-level pattern:

- Run variant calling
- Apply caller-specific filtering
- Annotate variants with VEP to generate a MAF

> **Caution:**
> Although scripts across callers are structurally similar, **inputs and
> assumptions differ slightly**. Always verify required arguments before running.

---

## 4.1.1 VarScan (Tumor-Only)

VarScan requires **two stages** due to how read count information is represented.

### Key Notes

- Tumor-only calling using `samtools mpileup`
- Read counts are **not reliably propagated** through VEP annotation
- Counts must be extracted and merged manually in Stage 2
- Uses the target BED file from **01 – Preparation**
- Requires:
  - `cnvkit_env` (utilities / merging)
  - `vep_env` (annotation)

---

### Stage 1 – Variant Calling, Filtering, and Annotation

This stage:
- Generates VarScan SNV calls
- Filters variants using a Panel of Normals (PoN)
- Annotates variants using VEP
- Produces an initial MAF per sample

#### Command

```bash
bash run_varscan_stage1.sh \
  <gatk_bam_config.yaml> \
  <pon_positions_file> \
  <targets.bed> \
  <reference.fa> \
  <output_directory>
```
---

### Stage 2 – Extract Read Counts and Finalize MAF

This stage:
- Extracts depth, reference, and alternate allele counts from the VCF
- Merges these counts into the annotated MAF
- Produces the **final VarScan MAF** used downstream

#### Command

```bash
bash run_varscan_stage2_counts.sh \
  <gatk_bam_config.yaml> \
  <output_directory>
```
---

### VarScan Outputs

For each sample, the final output is:

```text
<SAMPLE>_VarScan_annotated_with_counts.maf
```
This file is used in downstream multi-caller merging.

---

