#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_mutect2_annotate.sh bam_config.yaml ref.fa outdir

YAML="$1"
REF="$2"
OUTBASE="$3"

# Extract BAMs
BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)
    SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
    LOGDIR="$OUTBASE/$SAMPLE/logs"
    TEMP="$OUTBASE/$SAMPLE/TEMP"

    FILTERED="${SAMPLE_DIR}/${SAMPLE}_MuTect2_filtered.vcf"
    PASS="${SAMPLE_DIR}/${SAMPLE}_MuTect2_PASS.vcf"
    MAF="${SAMPLE_DIR}/${SAMPLE}_MuTect2_filtered_annotated.maf"
    ANNOTVCF="${SAMPLE_DIR}/${SAMPLE}_MuTect2_filtered_annotated.vcf"

    mkdir -p "$LOGDIR" "$TEMP"

    echo "Submitting Mutect2 annotation for $SAMPLE"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=mutect2_annot_${SAMPLE}
#SBATCH --output=${LOGDIR}/%x.out
#SBATCH --error=${LOGDIR}/%x.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH -c 4

set -euo pipefail

SAMPLE="${SAMPLE}"
SAMPLE_DIR="${SAMPLE_DIR}"
TEMP="${TEMP}"

FILTERED="\${SAMPLE_DIR}/\${SAMPLE}_MuTect2_filtered.vcf"
PASS="\${SAMPLE_DIR}/\${SAMPLE}_MuTect2_PASS.vcf"
MAF="\${SAMPLE_DIR}/\${SAMPLE}_MuTect2_filtered_annotated.maf"
ANNOTVCF="\${SAMPLE_DIR}/\${SAMPLE}_MuTect2_filtered_annotated.vcf"

module load vcftools/0.1.15
module load samtools/1.20
module load tabix

source /cluster/home/t138377uhn/miniconda/etc/profile.d/conda.sh
conda activate vep_env

echo "=== Filtering PASS variants ==="

vcftools \
    --vcf "\$FILTERED" \
    --stdout \
    --recode \
    --keep-filtered PASS \
    --temp "\$TEMP" \
    > "\$PASS"

md5sum "\$PASS" > "\${PASS}.md5"

echo "=== Running VEP + vcf2maf ==="

vcf2maf.pl \
  --species homo_sapiens \
  --ncbi-build GRCh38 \
  --ref-fasta "$REF" \
  --input-vcf "\$PASS" \
  --output-maf "\$MAF" \
  --tumor-id "\$SAMPLE" \
  --vep-path \$(dirname \$(which vep)) \
  --vep-data /cluster/projects/kridelgroup/resources/VEP/GRCh38/98 \
  --vep-forks 4 \
  --buffer-size 1000 \
  --tmp-dir "\$TEMP"

if [ -s "\$MAF" ]; then
    mv "\$TEMP/\${SAMPLE}_MuTect2_PASS.vep.vcf" "\$ANNOTVCF"
    md5sum "\$ANNOTVCF" > "\${ANNOTVCF}.md5"

    bgzip -f "\$ANNOTVCF"
    tabix -p vcf "\${ANNOTVCF}.gz"

    rm -rf "\$TEMP"
fi

echo "DONE Mutect2 annotate for \$SAMPLE"
EOF

done

echo "=== Submitted all Mutect2 annotation jobs ==="
