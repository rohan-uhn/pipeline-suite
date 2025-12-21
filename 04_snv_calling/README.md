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

## 4.1.1 VarScan

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

## 4.1.2 VarDict

VarDict is run as a **two-stage tumor-only variant calling workflow**, similar in
structure to VarScan, but with different internal filtering logic and output
characteristics.

---

### Key Notes

- Tumor-only calling using **VarDictJava**
- Uses the **same target BED file** generated in **01 – Preparation**
- Applies internal VarDict filters prior to PoN filtering
- Variants are filtered to SNVs and short indels
- Panel of Normals (PoN) is applied to remove recurrent artifacts
- Requires:
  - System-installed VarDictJava + R
  - `vep_env` for VEP annotation

---

## Stage 1 – Variant Calling and Filtering

This stage performs:
- Raw variant calling using VarDict
- Strand-bias correction and VCF conversion
- Filtering based on:
  - Variant depth
  - Strand bias
  - Variant type (SNVs and indels)
- Panel of Normals (PoN) exclusion
- Produces a filtered VCF per sample

### Command

```bash
bash run_vardict_stage1.sh \
  <gatk_bam_config.yaml> \
  <pon_positions_file> \
  <targets.bed> \
  <reference.fa> \
  <output_directory>
````

### Notes on Filtering

The following filtering logic is applied internally:
- Minimum variant allele frequency (default: 1%)
- Minimum supporting reads
- Strand bias filtering
- Variant type restricted to:
  - SNVs
  - Insertions
  - Deletions
- PoN filtering applied **after** VarDict-specific filters

The output of this stage is a **filtered VCF**, not yet annotated.

---

## Stage 2 – VEP Annotation and Final MAF Generation

This stage:
- Annotates the filtered VarDict VCF using VEP
- Converts annotated variants to MAF format using VEP

### Command

```bash
bash run_vardict_stage2.sh \
  <gatk_bam_config.yaml> \
  <reference.fa> \
  <output_directory>
```

---

## VarDict Outputs

For each sample, the final output is:

```text
<SAMPLE>_VarDict_filtered_annotated.maf
```
This file is used in downstream multi-caller merging.

---

## 4.1.3 Strelka

Strelka is run as a **two-stage tumor-only variant calling workflow** using
the **Strelka germline model**, restricted to targeted regions.

Although Strelka is primarily designed for germline or paired tumor–normal
analysis, it can be adapted for tumor-only SNV and indel calling by:
- Restricting calls to target regions
- Applying strict post-calling filtering
- Removing recurrent artifacts using a Panel of Normals (PoN)

---

### Key Notes
- Tumor-only calling using **Strelka germline workflow**
- Variant calling restricted to target regions
- Requires a **compressed and indexed BED file** (`.bed.gz` + `.tbi`)
  generated in **01 – Preparation**
- Requires a **compressed** Panel of Normals (PoN) (`.vcf.gz`)
- Requires:
  - Strelka
  - `vep_env` for annotation

---

## Stage 1 – Variant Calling and Filtering

This stage performs:
1. **Strelka germline variant calling**
   - Run on tumor-only BAMs
   - Restricted to targeted regions using the compressed BED file
2. **Post-calling filtering**
   - Retains only `PASS` variants
   - Removes variants present in the Panel of Normals (PoN)

The output of this stage is a **filtered Strelka VCF per sample**.

---

### Command

```bash
bash run_strelka_stage1.sh \
  <gatk_bam_config.yaml> \
  <panel_of_normals.vcf.gz> \
  <targets.bed.gz> \
  <reference.fa> \
  <output_directory>
````

> **Important:**
> The BED file and PoN VCF **must** be bgzipped and indexed (`.bed.gz` + `.tbi` + `vcf.gz`).
> The BED file is created during **01 – Preparation**.

---

## Stage 2 – VEP Annotation and Final MAF Generation

This stage:
- Annotates the filtered Strelka VCF using VEP
- Converts variants to MAF format using VEP
- Optionally retains an annotated VCF for reference and debugging
- Cleans up intermediate Strelka working directories after success

---

### Command

```bash
bash run_strelka_stage2_annotate.sh \
  <gatk_bam_config.yaml> \
  <reference.fa> \
  <output_directory>
```

Stage 2 must be run **after Stage 1 completes successfully**.

---

## Strelka Outputs

For each sample, the final output is:

```text
<SAMPLE>_Strelka_annotated.maf
```

---

## 4.1.4 MuTect

MuTect (v1) is run as a **two-stage tumor-only somatic SNV calling workflow**.
Although MuTect was originally designed for paired tumor–normal analysis, it
supports tumor-only calling through the use of a **Panel of Normals (PoN)** and
external variant databases.

---

### Key Notes
- Tumor-only somatic SNV calling using **MuTect v1**
- Requires a **MuTect-compatible Panel of Normals**
  - This PoN format differs from other callers
- Uses external variant resources:
  - **dbSNP (GRCh38 common variants)**
  - **COSMIC variants**
- Variant calling restricted to target regions
- Filtering removes MuTect `REJECT` calls prior to annotation
- Requires:
  - MuTect v1
  - `vep_env` for annotation

---

## Panel of Normals and Resource Files

MuTect1 requires additional reference resources beyond a standard PoN:
- **Panel of Normals (PoN)**  
  Used to remove recurrent technical artifacts in tumor-only mode
- **dbSNP (GRCh38)**  
  Used to annotate known common variants
- **COSMIC variants**  
  Used to annotate known somatic mutations

> **Important:**  
> For MuTect v1, minor VCF header modifications may be required:
> - Change the VCF version header from `VCFv4.2` to `VCFv4.1`
> - Change the `AD` INFO definition from `Number=R` to `Number=.`
>  
> These changes ensure compatibility with MuTect v1’s VCF parser.

---

## Stage 1 – Variant Calling and Filtering

This stage performs:
1. **MuTect tumor-only variant calling**
   - Runs MuTect with a tumor BAM only
   - Uses PoN, dbSNP, and COSMIC for artifact suppression and annotation
   - Restricts calling to target regions with padding

2. **Post-calling filtering**
   - Removes variants flagged as `REJECT`
   - Retains only high-confidence MuTect calls

The output of this stage is a **filtered MuTect VCF per sample**.

---

### Command

```bash
bash run_mutect_batch.sh \
  <gatk_bam_config.yaml> \
  <mutect_pon.vcf> \
  <dbsnp.vcf.gz> \
  <cosmic.vcf.gz> \
  <targets.bed> \
  <reference.fa> \
  <output_directory>
````

---

## Stage 2 – VEP Annotation and Final MAF Generation

This stage:
- Annotates the filtered MuTect VCF using VEP
- Converts variants to MAF format using VEP
- Optionally retains an annotated VCF for reference
- Cleans up intermediate files upon successful annotation

---

### Command

```bash
bash run_mutect_annotate.sh \
  <gatk_bam_config.yaml> \
  <reference.fa> \
  <output_directory>
```

Stage 2 must be run **after Stage 1 completes successfully**.

---

## MuTect Outputs

For each sample, the final output is:

```text
<SAMPLE>_MuTect_filtered_annotated.maf
```

---

