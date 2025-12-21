import pandas as pd
import os

# List of sample IDs
with open("list.txt") as f:
    sample_ids = [line.strip() for line in f if line.strip()]

# Shared columns (excluding the ones to drop)
shared_columns = [
    "Hugo_Symbol", "Entrez_Gene_Id", "NCBI_Build", "Chromosome", "Start_Position", "End_Position",
    "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "all_effects", "Allele",
    "Gene", "Feature", "Feature_type", "Consequence", "Existing_variation"
]

# Caller-specific numeric fields
caller_specific_cols = ["t_depth", "t_ref_count", "t_alt_count"]

# Unique variant identifier
id_columns = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]

# Callers to loop through
callers = ["MuTect2", "MuTect", "VarDict", "VarScan", "Strelka"]

for sample in sample_ids:
    print(f"\n Processing sample: {sample}")

    # Build file paths
    caller_paths = {
        "MuTect2": f"/cluster/projects/kridelgroup/Rohan/finalized_calling/TROG/Mutect2/{sample}/{sample}/{sample}_MuTect2_filtered_annotated.maf",
        "MuTect": f"/cluster/projects/kridelgroup/Rohan/finalized_calling/TROG/Mutect/{sample}/{sample}/{sample}_MuTect_filtered_annotated.maf",
        "VarDict": f"/cluster/projects/kridelgroup/Rohan/finalized_calling/TROG/Vardict/{sample}/{sample}/{sample}_VarDict_filtered_annotated.maf",
        "VarScan": f"/cluster/projects/kridelgroup/Rohan/finalized_calling/TROG/Varscan/{sample}/{sample}/{sample}_VarScan_annotated_with_counts.maf",
        "Strelka": f"//cluster/projects/kridelgroup/Rohan/finalized_calling/TROG/Strelka/{sample}/{sample}/Strelka/{sample}_Strelka_annotated.maf"
    }

    all_dfs = []
    for caller, path in caller_paths.items():
        if os.path.exists(path):
            print(f"Loading {caller} for {sample}")
            df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
            df = df.loc[:, ~df.columns.duplicated()]  # Remove duplicates

            # Ensure required columns are present
            for col in shared_columns + caller_specific_cols:
                if col not in df.columns:
                    df[col] = None

            df = df[shared_columns + caller_specific_cols + id_columns].copy()

            # Rename caller-specific fields and drop originals
            for col in caller_specific_cols:
                df[f"{col}_{caller}"] = df[col]
                df.drop(columns=[col], inplace=True)

            df[caller] = True
            df = df.loc[:, ~df.columns.duplicated()]
            df = df.drop_duplicates(subset=id_columns)
            all_dfs.append(df)
        else:
            print(f"File not found: {path}")

    if not all_dfs:
        print(f"No files found for {sample}, skipping.")
        continue

    # Combine all variant calls from different callers
    merged_df = pd.concat(all_dfs, axis=0, ignore_index=True)

    # Group by unique variant and aggregate shared columns
    final_df = merged_df.groupby(id_columns, dropna=False).agg(
        lambda x: x.dropna().iloc[0] if x.notna().any() else None
    ).reset_index()


    # Parse Existing_variation into rsIDs, COSMIC IDs, and other IDs
    def split_existing_variation(ev):
        if pd.isna(ev):
            return pd.Series([None, None, None])
        rs = []
        cosm = []
        other = []
        for entry in str(ev).split(","):
            entry = entry.strip()
            if entry.startswith("rs"):
                rs.append(entry)
            elif entry.startswith("COSM"):
                cosm.append(entry)
            elif entry:
                other.append(entry)
        return pd.Series([
            ",".join(rs) if rs else None,
            ",".join(cosm) if cosm else None,
            ",".join(other) if other else None
        ])

    # Apply the splitting logic
    final_df[["Existing_rsIDs", "Existing_COSMIC", "Existing_other"]] = final_df["Existing_variation"].apply(split_existing_variation)


    # Ensure caller presence flags are consistent
    for caller in callers:
        if caller not in final_df.columns:
            final_df[caller] = False
        else:
            final_df[caller] = merged_df.groupby(id_columns, dropna=False)[caller].max().reset_index(drop=True).fillna(False)

    # Ensure caller-specific columns exist across all callers
    for caller in callers:
        for col in caller_specific_cols:
            colname = f"{col}_{caller}"
            if colname not in final_df.columns:
                final_df[colname] = None

    # Save output
    output_path = f"merged_{sample}_variants_cleaned.tsv"
    final_df.to_csv(output_path, sep="\t", index=False)
    print(f"Saved: {output_path}")
