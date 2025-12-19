set -euo pipefail

# CHECK ARGUMENTS
if [ "$#" -ne 3 ]; then
    echo "Usage:"
    echo "  sbatch run_tp53_cnv_wrapper.sh <sample_list.txt> <bamqc_dir> <outdir>"
    exit 1
fi

SAMPLES="$1"
BAMQC="$2"
OUTDIR="$3"

echo "=== TP53 CNV WRAPPER ==="
echo "SAMPLES:  $SAMPLES"
echo "BAMQC:    $BAMQC"
echo "OUTDIR:   $OUTDIR"
echo "========================="
echo

# LOAD R
module load R/3.6.1 2>/dev/null || module load R 2>/dev/null

# RUN THE PIPELINE
echo "Running TP53 CNV pipeline..."

Rscript tp53_cnv_from_bamqc_hg38_from_list.R \
    "$SAMPLES" \
    "$BAMQC" \
    "$OUTDIR"

echo
echo "=== TP53 CNV pipeline finished successfully ==="
echo "Results saved to: $OUTDIR"
