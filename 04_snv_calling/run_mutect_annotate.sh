#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_mutect_annotate.sh bam_config.yaml ref.fa outdir

YAML="$1"
REF="$2"
OUTBASE="$3"

# Extract BAMs
BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)

    SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
    LOGDIR="$OUTBASE/$SAMPLE/logs"
    TEMP_DIR="$SAMPLE_DIR/TEMP"

    FILTERED_VCF="$SAMPLE_DIR/${SAMPLE}_MuTect_filtered.vcf"
    MAF="$SAMPLE_DIR/${SAMPLE}_MuTect_filtered_annotated.maf"
    ANNOT_VCF="$SAMPLE_DIR/${SAMPLE}_MuTect_filtered_annotated.vcf"

    mkdir -p "$LOGDIR"

    echo "Submitting Mutect1 annotation for $SAMPLE"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=annot_mutect_${SAMPLE}
#SBATCH --output=${LOGDIR}/%x.out
#SBATCH --error=${LOGDIR}/%x.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH -c 4

set -euo pipefail

module load samtools/1.20
module load tabix

source /cluster/home/t138377uhn/miniconda/etc/profile.d/conda.sh
conda activate vep_env

echo "=== Annotating MuTect1 for $SAMPLE ==="

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
    echo "[SUCCESS] MAF exists, finalizing annotation for $SAMPLE"

    if [[ -s "$TEMP_DIR/${SAMPLE}_MuTect_filtered.vep.vcf" ]]; then
        mv "$TEMP_DIR/${SAMPLE}_MuTect_filtered.vep.vcf" "$ANNOT_VCF"
        md5sum "$ANNOT_VCF" > "$ANNOT_VCF.md5"

        bgzip -f "$ANNOT_VCF"
        tabix -p vcf "$ANNOT_VCF.gz"
    fi

    rm -rf "$TEMP_DIR"
else
    echo "ERROR: MAF missing for $SAMPLE — not cleaning TEMP"
fi

echo "=== DONE annotation for $SAMPLE ==="
EOF

done

echo "=== Submitted all Mutect annotation jobs ==="
