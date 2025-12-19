#!/bin/bash
set -euo pipefail

# Usage:
#   bash make_cnvkit_reference_cleaned.sh bam_config.yaml raw_targets.bed outdir genome.fa

YAML="$1"     
BED="$2"
OUTDIR="$3"
REF="$4"

mkdir -p "$OUTDIR"

echo "=== CNVkit Reference Builder (with BED cleaning) ==="
echo "BED input:      $BED"
echo "Output dir:     $OUTDIR"
echo "Reference FASTA: $REF"
echo

# Clean BED file
module load bedtools

FAI="${REF}.fai"
CLEAN_BED="$OUTDIR/clean_targets.bed"

echo "[1/3] Cleaning BED file..."
echo "    Output: $CLEAN_BED"

awk 'BEGIN{OFS="\t"} 
    $0!~/^track|^#/ && NF>=3 && $2<$3 {print $0}' "$BED" \
| bedtools sort -faidx "$FAI" -i - \
| awk 'BEGIN{OFS="\t"} {key=$1 FS $2 FS $3} !seen[key]++' \
> "$CLEAN_BED"

echo "    BED cleaned and sorted."
echo

# Create target bins
TARGETS="$OUTDIR/targets.target.bed"

echo "[2/3] Creating target bins..."
cnvkit.py target "$CLEAN_BED" \
    -o "$TARGETS"

echo "    Targets saved to $TARGETS"
echo

# Build flat CNVkit reference (target-only)
REFERENCE="$OUTDIR/flat_reference.cnn"

echo "[3/3] Building flat reference..."
cnvkit.py reference \
    -t "$TARGETS" \
    -f "$REF" \
    -o "$REFERENCE"

echo "=== Done. CNVkit reference created ==="
echo "Reference: $REFERENCE"
echo "Clean BED: $CLEAN_BED"
echo "Targets:   $TARGETS"
