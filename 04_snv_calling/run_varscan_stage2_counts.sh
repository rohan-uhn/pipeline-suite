#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_varscan_stage2_counts.sh bam_config.yaml outdir

YAML="$1"
OUTBASE="$2"

# Extract sample names from YAML
BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)
    SAMPLE_DIR="$OUTBASE/$SAMPLE/$SAMPLE"
    LOGDIR="$OUTBASE/$SAMPLE/logs_stage2"
    mkdir -p "$LOGDIR"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=varscan2_$SAMPLE
#SBATCH --output=$LOGDIR/%x.out
#SBATCH --error=$LOGDIR/%x.err
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -c 2

set -euo pipefail
echo "Running Stage 2 VarScan counts merge for $SAMPLE"

module load samtools/1.20
module load vcftools/0.1.15
module load tabix

# Activate clean Python environment
source /cluster/home/t138377uhn/miniconda/etc/profile.d/conda.sh
conda activate cnvkit_env

FILTERED_VCF="$SAMPLE_DIR/${SAMPLE}_VarScan_filtered.vcf"
COUNTS="$SAMPLE_DIR/vcf_counts.tsv"
MAF="$SAMPLE_DIR/${SAMPLE}_VarScan_annotated.maf"
OUT_MAF="$SAMPLE_DIR/${SAMPLE}_VarScan_annotated_with_counts.maf"

# 1) EXTRACT DP / REF / ALT COUNTS

if bcftools view -h "\$FILTERED_VCF" | grep -q 'ID=RD'; then
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%RD]\t[%AD]\n' "\$FILTERED_VCF" > "\$COUNTS"

elif bcftools view -h "\$FILTERED_VCF" | grep -q 'ID=AD'; then
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\n' "\$FILTERED_VCF" \
      | awk -F'\t' 'BEGIN{OFS="\t"}{
            split(\$6,a,","); ref=a[1]; alt=a[2];
            d=\$5; if(d=="."||d=="") d = ref+alt;
            print \$1,\$2,\$3,\$4,d,ref,alt
        }' > "\$COUNTS"

elif bcftools view -h "\$FILTERED_VCF" | grep -q 'ID=DP4'; then
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP4]\n' "\$FILTERED_VCF" \
      | awk -F'\t' 'BEGIN{OFS="\t"}{
            split(\$5,a,","); ref=a[1]+a[2]; alt=a[3]+a[4];
            print \$1,\$2,\$3,\$4,ref+alt,ref,alt
      }' > "\$COUNTS"

else
    echo "WARN: No count fields. Writing placeholders."
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t.\t.\t.\n' "\$FILTERED_VCF" > "\$COUNTS"
fi

# 2) MERGE COUNTS INTO MAF (Python)

python3 <<PYCODE
import pandas as pd
import os

sample = "${SAMPLE}"
sample_dir = "${SAMPLE_DIR}"

maf_path = os.path.join(sample_dir, sample + "_VarScan_annotated.maf")
vcf_path = os.path.join(sample_dir, "vcf_counts.tsv")
out_path = os.path.join(sample_dir, sample + "_VarScan_annotated_with_counts.maf")

maf = pd.read_csv(maf_path, sep="\t", comment="#", low_memory=False)
vcf = pd.read_csv(vcf_path, sep="\t", header=None)
vcf.columns=["Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2",
             "t_depth","t_ref_count","t_alt_count"]

maf["Start_Position"] = maf["Start_Position"].astype(int)
vcf["Start_Position"] = vcf["Start_Position"].astype(int)

merged = pd.merge(
    maf, vcf,
    on=["Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2"],
    how="left"
)

merged.to_csv(out_path, sep="\t", index=False)
print("WROTE:", out_path)
PYCODE

echo "Finished Stage 2 for $SAMPLE"
EOF

done
