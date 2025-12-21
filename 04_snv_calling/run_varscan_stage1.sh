#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_varscan_stage1.sh bam_config.yaml pon.vcf bedfile.bed ref.fa outdir

YAML="$1"
PON="$2"
BED="$3"
REF="$4"
OUTBASE="$5"

mkdir -p "$OUTBASE"

# Extract BAM paths from YAML
BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)
    SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
    LOGDIR="$OUTBASE/$SAMPLE/logs_stage1"

    mkdir -p "$SAMPLE_DIR" "$LOGDIR"

    echo "Submitting VarScan Stage 1 pipeline for $SAMPLE"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=varscan1_$SAMPLE
#SBATCH --output=$LOGDIR/%x.out
#SBATCH --error=$LOGDIR/%x.err
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -c 4

set -euo pipefail
echo "Running Stage 1 VarScan for $SAMPLE"

module load samtools/1.20
module load varscan/2.4.2
module load vcftools/0.1.15
module load tabix

# Activate VEP environment for annotation
source /cluster/home/t138377uhn/miniconda/etc/profile.d/conda.sh
conda activate vep_env

RAW_VCF="$SAMPLE_DIR/${SAMPLE}_VarScan_raw.vcf"
FILTERED_VCF="$SAMPLE_DIR/${SAMPLE}_VarScan_filtered.vcf"
MAF="$SAMPLE_DIR/${SAMPLE}_VarScan_annotated.maf"
TEMP_DIR="$SAMPLE_DIR/TEMP"

mkdir -p "\$TEMP_DIR"

# 1) RUN VARSCAN

samtools mpileup -B -q1 -d10000 -f "$REF" -l "$BED" "$BAM" \
  | awk -F"\t" '\$4 > 0' \
  | java -Xmx10g -Djava.io.tmpdir=\$TEMP_DIR \
      -jar \$varscan_dir/VarScan.jar mpileup2cns - \
      --output-vcf 1 --variants 1 \
      > "\$RAW_VCF"

md5sum "\$RAW_VCF" > "\$RAW_VCF.md5"

# 2) FILTER USING PANEL OF NORMALS

vcftools --vcf "\$RAW_VCF" \
  --keep-filtered PASS \
  --stdout \
  --recode \
  --temp "\$TEMP_DIR" \
  --exclude-positions "$PON" \
  > "\$FILTERED_VCF"

md5sum "\$FILTERED_VCF" > "\$FILTERED_VCF.md5"

# 3) RUN VEP + VCF2MAF

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

echo "Finished Stage 1 for $SAMPLE"
EOF

done
