#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_mutect2.sh bam_config.yaml pon.vcf bedfile.bed ref.fa outdir

YAML="$1"
PON="$2"
BED="$3"
REF="$4"
OUTBASE="$5"

mkdir -p "$OUTBASE"

# Extract BAMs
BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)
    SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
    LOGDIR="$OUTBASE/$SAMPLE/logs"
    TEMP="$OUTBASE/$SAMPLE/TEMP"

    mkdir -p "$SAMPLE_DIR" "$TEMP" "$LOGDIR"

    echo "Submitting Mutect2 job for $SAMPLE"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=mutect2_${SAMPLE}
#SBATCH --output=${LOGDIR}/%x.out
#SBATCH --error=${LOGDIR}/%x.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH -c 2

set -euo pipefail

module load gatk/4.6.0.0

RAW_VCF="${SAMPLE_DIR}/${SAMPLE}_MuTect2_raw.vcf"
FILTERED_VCF="${SAMPLE_DIR}/${SAMPLE}_MuTect2_filtered.vcf"

echo "=== Running Mutect2 for $SAMPLE ==="

gatk Mutect2 \
    -R "$REF" \
    -I "$BAM" \
    --panel-of-normals "$PON" \
    -O "\$RAW_VCF" \
    --intervals "$BED" \
    --interval-padding 100

md5sum "\$RAW_VCF" > "\${RAW_VCF}.md5"

echo "=== Filtering Mutect Calls ==="

gatk FilterMutectCalls \
    -R "$REF" \
    -V "\$RAW_VCF" \
    -O "\$FILTERED_VCF"

md5sum "\$FILTERED_VCF" > "\${FILTERED_VCF}.md5"

echo "DONE Mutect2 for $SAMPLE"
EOF

done

echo "=== Submitted all Mutect2 jobs ==="
