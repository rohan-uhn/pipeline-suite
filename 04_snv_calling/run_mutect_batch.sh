#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_mutect_batch.sh bam_config.yaml pon.vcf dbsnp.vcf cosmic.vcf bedfile.bed ref.fa outdir

YAML="$1"
PON="$2"
DBSNP="$3"
COSMIC="$4"
BED="$5"
REF="$6"
OUTBASE="$7"

mkdir -p "$OUTBASE"

# Extract BAM paths
BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)

    SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
    LOGDIR="$OUTBASE/$SAMPLE/logs"

    mkdir -p "$SAMPLE_DIR" "$LOGDIR"

    RAW_VCF="$SAMPLE_DIR/${SAMPLE}_MuTect.vcf"
    STATS="$SAMPLE_DIR/${SAMPLE}_MuTect.stats"
    FILTERED_VCF="$SAMPLE_DIR/${SAMPLE}_MuTect_filtered.vcf"
    TEMP_DIR="$SAMPLE_DIR/TEMP"

    mkdir -p "$TEMP_DIR"

    echo "Submitting MuTect1 job for $SAMPLE"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=mutect_${SAMPLE}
#SBATCH --output=${LOGDIR}/%x.out
#SBATCH --error=${LOGDIR}/%x.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH -c 2

set -euo pipefail

module load mutect/1.1.5
module load vcftools/0.1.15

echo "=== Running MuTect1 for $SAMPLE ==="

java -Xmx12g -Djava.io.tmpdir=$TEMP_DIR \
  -jar \$mutect_dir/muTect.jar \
  -T MuTect \
  -R "$REF" \
  --input_file:tumor "$BAM" \
  --tumor_sample_name "$SAMPLE" \
  --vcf "$RAW_VCF" \
  --out "$STATS" \
  --normal_panel "$PON" \
  --dbsnp "$DBSNP" \
  --cosmic "$COSMIC" \
  --intervals "$BED" \
  --interval_padding 100

md5sum "$RAW_VCF" > "$RAW_VCF.md5"

echo "=== Filtering MuTect output for $SAMPLE ==="

vcftools \
  --vcf "$RAW_VCF" \
  --remove-filtered REJECT \
  --stdout \
  --recode \
  --temp "$TEMP_DIR" \
  > "$FILTERED_VCF"

md5sum "$FILTERED_VCF" > "$FILTERED_VCF.md5"

echo "Done MuTect1 + filtering for $SAMPLE"
EOF

done

echo "=== Submitted all MuTect jobs ==="
