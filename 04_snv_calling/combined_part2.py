import pandas as pd

# Path to list of files
with open("caller_list.txt") as f:
    file_list = [line.strip() for line in f if line.strip()]

# Output file name
output_path = "merged_all_samples_variants_cleaned.tsv"

# Initialize an empty list to store DataFrames
combined_df = None

for i, file in enumerate(file_list):
    print(f"Adding: {file}")
    df = pd.read_csv(file, sep="\t", low_memory=False)

    if combined_df is None:
        combined_df = df  # First file includes header
    else:
        # Subsequent files: make sure columns match and only append rows
        df = df[combined_df.columns.tolist()]
        combined_df = pd.concat([combined_df, df], axis=0, ignore_index=True)

# Save the final combined file
combined_df.to_csv(output_path, sep="\t", index=False)
print(f"\nCombined file saved to: {output_path}")
