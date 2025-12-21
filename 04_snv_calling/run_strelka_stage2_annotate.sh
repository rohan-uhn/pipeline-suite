#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_strelka_stage2_annotate.sh \
#     bam_config.yaml \
#     reference.fa \
#     outdir

YAML="$1"
REF="$2"
OUTBASE="$3"

# Extract BAMs 
mapfile -t BAMS < <(
  awk '
    /recalibrated\.bam$/ {
      for (i=1; i<=NF; i++) {
        if ($i ~ /recalibrated\.bam$/) print $i
      }
    }
  ' "$YAML"
)

if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "ERROR: No BAMs found in YAML" >&2
  exit 1
fi

# Loop over samples
for BAM in "${BAMS[@]}"; do

  SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)

  SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
  STRELKA_DIR="$SAMPLE_DIR/Strelka"
  TEMP_DIR="$SAMPLE_DIR/TEMP"
  LOGDIR="$OUTBASE/$SAMPLE/logs_stage2"

  FILTERED_VCF="$STRELKA_DIR/${SAMPLE}_Strelka_filtered.vcf"
  MAF="$STRELKA_DIR/${SAMPLE}_Strelka_annotated.maf"
  ANNOT_VCF="$STRELKA_DIR/${SAMPLE}_Strelka_annotated.vcf"

  mkdir -p "$LOGDIR"

  echo "Submitting Strelka annotation job for $SAMPLE"

  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=annot_strelka_${SAMPLE}
#SBATCH --output=${LOGDIR}/%x.out
#SBATCH --error=${LOGDIR}/%x.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH -c 4

set -euo pipefail

echo "=== Annotating Strelka for $SAMPLE ==="
echo "HOSTNAME: \$HOSTNAME"

module load samtools/1.20
module load tabix

source /cluster/home/t138377uhn/miniconda/etc/profile.d/conda.sh
conda activate vep_env

if [[ ! -s "$FILTERED_VCF" ]]; then
  echo "ERROR: Filtered VCF not found for $SAMPLE"
  exit 1
fi

vcf2maf.pl \
  --species homo_sapiens \
  --ncbi-build GRCh38 \
  --ref-fasta "$REF" \
  --input-vcf "$FILTERED_VCF" \
  --output-maf "$MAF" \
  --tumor-id "$SAMPLE" \
  --vep-path \$(dirname \$(which vep)) \
  --vep-data /cluster/projects/kridelgroup/resources/VEP/GRCh38/98 \
  --vep-forks 4 \
  --buffer-size 1000 \
  --tmp-dir "$TEMP_DIR"

if [[ -s "$MAF" ]]; then
  echo "[SUCCESS] MAF generated for $SAMPLE"

  if [[ -s "$TEMP_DIR/${SAMPLE}_Strelka_filtered.vep.vcf" ]]; then
    mv "$TEMP_DIR/${SAMPLE}_Strelka_filtered.vep.vcf" "$ANNOT_VCF"
    md5sum "$ANNOT_VCF" > "$ANNOT_VCF.md5"

    bgzip -f "$ANNOT_VCF"
    tabix -p vcf "$ANNOT_VCF.gz"
  else
    echo "WARNING: Annotated VCF not found in TEMP"
  fi

  rm -rf "$TEMP_DIR"
  rm -rf "$STRELKA_DIR/workspace" "$STRELKA_DIR/results"

else
  echo "ERROR: MAF missing for $SAMPLE — skipping cleanup."
fi

echo "=== DONE annotation for $SAMPLE ==="
EOF

done

echo "=== Submitted all Strelka Stage 2 annotation jobs ==="
