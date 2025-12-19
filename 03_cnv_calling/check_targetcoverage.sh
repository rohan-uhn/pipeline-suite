#!/bin/bash
set -euo pipefail

# Usage:
#   bash check_targetcoverage.sh /path/to/outdir

OUTDIR="$1"

echo "Checking for *.cnr in each sample directory..."
echo

missing=0
total=0

for SAMPLE_DIR in "$OUTDIR"/*/; do
    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")
    total=$((total + 1))

    # Look for ANY file ending with .cnr
    COVER=$(ls "$SAMPLE_DIR"/*.cnr 2>/dev/null | wc -l)

    if [[ "$COVER" -gt 0 ]]; then
        echo "[OK]   $SAMPLE"
    else
        echo "[MISS] $SAMPLE --> No .cnr file found"
        missing=$((missing + 1))
    fi
done

echo
echo "Summary: $missing missing out of $total sample directories."
