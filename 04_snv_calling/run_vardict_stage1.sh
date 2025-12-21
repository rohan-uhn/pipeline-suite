#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_vardict_stage1.sh bam_config.yaml pon.vcf bedfile.bed ref.fa outdir

YAML="$1"
PON="$2"
BED="$3"
REF="$4"
OUTBASE="$5"

mkdir -p "$OUTBASE"

# Extract BAMs from YAML
BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)
    SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
    LOGDIR="$OUTBASE/$SAMPLE/logs_stage1"
    mkdir -p "$SAMPLE_DIR" "$LOGDIR"

    echo "Submitting VarDict Stage 1 for $SAMPLE"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=vardict1_$SAMPLE
#SBATCH --output=$LOGDIR/%x.out
#SBATCH --error=$LOGDIR/%x.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH -c 4

set -euo pipefail
echo "Running VarDict Stage 1 for $SAMPLE"
echo "Host: \$HOSTNAME"

module load perl
module load vardictjava/1.7.0
module load R/4.1.0
module load samtools/1.20
module load vcftools/0.1.15
module load tabix

DIRNAME=\$(dirname \$(which VarDict))
export JAVA_OPTS="-Xmx15g"

RAW_VCF="$SAMPLE_DIR/${SAMPLE}_VarDict_raw.vcf"
FILTERED_VCF="$SAMPLE_DIR/${SAMPLE}_VarDict_filtered.vcf"
TEMP_DIR="$SAMPLE_DIR/TEMP"

mkdir -p "\$TEMP_DIR"

# 1) RUN VARDICT

VarDict \
  -G "$REF" \
  -f 0.01 \
  -U \
  -N "$SAMPLE" \
  -c 1 -S 2 -E 3 \
  "$BED" \
  -b "$BAM" \
 | \$DIRNAME/teststrandbias.R \
 | \$DIRNAME/var2vcf_valid.pl -N "$SAMPLE" -S -f 0.01 \
 > "\$RAW_VCF"

md5sum "\$RAW_VCF" > "\$RAW_VCF.md5"

# 2) FILTER RAW VARDICT VCF

# Compress/index if needed
if [ ! -s "\$RAW_VCF.gz.tbi" ]; then
    bgzip -f "\$RAW_VCF"
    tabix -p vcf "\$RAW_VCF.gz"
fi

bcftools filter --include 'INFO/SBF>0.1 & INFO/VD>10' "\$RAW_VCF.gz" \
  | perl -e '
        while(<>) {
            if(/^#/) { print \$_; next; }
            @a = split /\\t/;
            if(\$a[7] =~ /TYPE=SNV|TYPE=Insertion|TYPE=Deletion/) {
                print join("\\t", @a);
            }
        }
    ' \
  | vcftools \
        --vcf - \
        --keep-filtered PASS \
        --stdout \
        --recode --recode-INFO-all \
        --exclude-positions "$PON" \
        --temp "\$TEMP_DIR" \
  > "\$FILTERED_VCF"

md5sum "\$FILTERED_VCF" > "\$FILTERED_VCF.md5"

echo "Stage 1 completed for $SAMPLE"
EOF

done
