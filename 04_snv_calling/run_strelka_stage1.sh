#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_strelka_stage1.sh \
#     bam_config.yaml \
#     panel_of_normals.vcf[.gz] \
#     targets.bed.gz \
#     reference.fa \
#     outdir
#

YAML="$1"
PON="$2"
BED="$3"
REF="$4"
OUTBASE="$5"

mkdir -p "$OUTBASE"

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

echo "Found ${#BAMS[@]} BAM(s)"

# Loop over samples
for BAM in "${BAMS[@]}"; do

  SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)

  SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
  STRELKA_DIR="$SAMPLE_DIR/Strelka"
  TEMP_DIR="$SAMPLE_DIR/TEMP"

  LOGDIR="$OUTBASE/$SAMPLE/logs_stage1"
  FILTER_LOG="$LOGDIR/filter"

  mkdir -p "$STRELKA_DIR" "$TEMP_DIR" "$LOGDIR" "$FILTER_LOG"

  echo "Submitting Strelka Stage 1 for $SAMPLE"

  # 1) STRELKA CALLING
  STRELKA_JOBID=$(sbatch <<EOF | awk '{print $4}'
#!/bin/bash
#SBATCH --job-name=run_strelka_${SAMPLE}
#SBATCH -D $LOGDIR
#SBATCH -t 02-00:00:00
#SBATCH --mem 8G
#SBATCH -c 1

set -e
echo HOSTNAME: \$HOSTNAME
echo SLURM NODE: \$SLURMD_NODENAME

module load strelka/2.9.10
module load python2/2.7.15

if [[ ! -s $STRELKA_DIR/runWorkflow.py ]]; then
  configureStrelkaGermlineWorkflow.py \
    --bam $BAM \
    --referenceFasta $REF \
    --runDir $STRELKA_DIR \
    --targeted \
    --callRegions $BED
fi

if [[ -s $STRELKA_DIR/workspace/pyflow.data/active_pyflow_process.txt ]]; then
  rm $STRELKA_DIR/workspace/pyflow.data/active_pyflow_process.txt
fi

$STRELKA_DIR/runWorkflow.py --quiet -m local
EOF
)

  # 2) FILTERING (vcftools logic)
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=filter_strelka_${SAMPLE}
#SBATCH -D $FILTER_LOG
#SBATCH -t 02-00:00:00
#SBATCH --mem 2G
#SBATCH -c 1
#SBATCH --dependency=afterok:${STRELKA_JOBID}
#SBATCH --kill-on-invalid-dep=yes

set -e
echo HOSTNAME: \$HOSTNAME
echo SLURM NODE: \$SLURMD_NODENAME

module load vcftools/0.1.15

vcftools \
  --gzvcf $STRELKA_DIR/results/variants/variants.vcf.gz \
  --stdout \
  --recode \
  --keep-filtered PASS \
  --temp $TEMP_DIR \
  --exclude-positions $PON \
  > $STRELKA_DIR/${SAMPLE}_Strelka_filtered.vcf

md5sum $STRELKA_DIR/${SAMPLE}_Strelka_filtered.vcf \
  > $STRELKA_DIR/${SAMPLE}_Strelka_filtered.vcf.md5
EOF

done

echo "=== Submitted all Strelka Stage 1 jobs ==="
