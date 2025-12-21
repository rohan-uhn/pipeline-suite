#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_vardict_stage2.sh bam_config.yaml ref.fa outdir

YAML="$1"
REF="$2"
OUTBASE="$3"

BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)
    SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
    LOGDIR="$OUTBASE/$SAMPLE/logs_stage2"
    mkdir -p "$LOGDIR"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=vardict2_$SAMPLE
#SBATCH --output=$LOGDIR/%x.out
#SBATCH --error=$LOGDIR/%x.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH -c 4

set -euo pipefail
echo "Running VarDict Stage 2 for $SAMPLE"

module load samtools/1.20
module load tabix

# Use ONLY VEP CONDA ENVIRONMENT
source /cluster/home/t138377uhn/miniconda/etc/profile.d/conda.sh
conda activate vep_env

FILTERED_VCF="$SAMPLE_DIR/${SAMPLE}_VarDict_filtered.vcf"
MAF="$SAMPLE_DIR/${SAMPLE}_VarDict_filtered_annotated.maf"
ANNOT_VCF="$SAMPLE_DIR/${SAMPLE}_VarDict_filtered_annotated.vcf"
TEMP_DIR="$SAMPLE_DIR/TEMP"

mkdir -p "\$TEMP_DIR"

# 1) RUN VEP + VCF2MAF

vcf2maf.pl \
  --species homo_sapiens \
  --ncbi-build GRCh38 \
  --ref-fasta "$REF" \
  --input-vcf "\$FILTERED_VCF" \
  --output-maf "\$MAF" \
  --tumor-id "$SAMPLE" \
  --vep-path \$(dirname \$(which vep)) \
  --vep-data /cluster/projects/kridelgroup/resources/VEP/GRCh38/98 \
  --vep-forks 4 \
  --buffer-size 1000 \
  --tmp-dir "\$TEMP_DIR"

# Move annotated VCF if present
if [ -s "\$TEMP_DIR/${SAMPLE}_VarDict_filtered.vep.vcf" ]; then
    mv "\$TEMP_DIR/${SAMPLE}_VarDict_filtered.vep.vcf" "\$ANNOT_VCF"
    md5sum "\$ANNOT_VCF" > "\$ANNOT_VCF.md5"

    bgzip -f "\$ANNOT_VCF"
    tabix -p vcf "\$ANNOT_VCF.gz"
fi

echo "Stage 2 completed for $SAMPLE"
EOF

done
