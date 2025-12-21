import pandas as pd
import pysam

# Load your combined TSV
df = pd.read_csv("merged_all_samples_variants_cleaned.tsv", sep="\t")

# Path to sorted and indexed gnomAD VCF
gnomad_vcf_path = "/cluster/projects/kridelgroup/Rohan/pugh_pipeline/transfered_references/hg38/af-only-gnomad.sorted.hg38.vcf.gz"
gnomad_vcf = pysam.TabixFile(gnomad_vcf_path)

# Function to get gnomAD AF from VCF for each row
def get_gnomad_af(chrom, pos, ref, alt):
    chrom = "chr" + chrom if not chrom.startswith("chr") else chrom
    try:
        records = gnomad_vcf.fetch(chrom, pos - 1, pos)
        for rec in records:
            fields = rec.strip().split('\t')
            vcf_ref = fields[3]
            vcf_alt = fields[4]
            info = fields[7]

            if vcf_ref == ref and vcf_alt == alt:
                for entry in info.split(";"):
                    if entry.startswith("AF="):
                        return float(entry.split("=")[1])
    except Exception:
        return None
    return None

# Ensure proper formatting of chromosome values
df["Chromosome"] = df["Chromosome"].astype(str).str.replace("^chr", "", regex=True)

# Apply function row-wise
df["gnomAD_AF"] = df.apply(
    lambda row: get_gnomad_af(row["Chromosome"], row["Start_Position"], row["Reference_Allele"], row["Tumor_Seq_Allele2"]),
    axis=1
)

# Save updated TSV
df.to_csv("merged_all_samples_with_gnomAD_AF.tsv", sep="\t", index=False)
print("gnomAD_AF column added.")
