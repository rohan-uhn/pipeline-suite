#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_cnvkit_per_sample.sh \
#       bam_config.yaml \
#       targets.target.bed \
#       flat_reference.cnn \
#       output_dir

YAML="$1"
TARGETS="$2"       
REFERENCE="$3"
OUTDIR="$4"

mkdir -p "$OUTDIR"

echo "=== [CNVkit] Running per-sample tumor-only CNVkit (TARGET-ONLY MODE) ==="
echo "YAML:        $YAML"
echo "TARGETS:     $TARGETS"
echo "REFERENCE:   $REFERENCE"
echo "OUTDIR:      $OUTDIR"
echo "==============================================================="
echo

# Extract BAM paths from YAML (all lines ending in recalibrated.bam)
BAMS=$(awk '/recalibrated\.bam$/ {print $2}' "$YAML")

for BAM in $BAMS; do
    SAMPLE=$(basename "$BAM" _realigned_recalibrated.bam)
    SAMPLE_DIR="$OUTDIR/$SAMPLE"
    LOGDIR="$SAMPLE_DIR/logs"

    mkdir -p "$SAMPLE_DIR" "$LOGDIR"

    echo "Submitting CNVkit job for: $SAMPLE"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=cnvkit_$SAMPLE
#SBATCH --output=$LOGDIR/%x.out
#SBATCH --error=$LOGDIR/%x.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

source ~/.bashrc
conda activate cnvkit_env

echo "Running CNVkit for sample $SAMPLE"

cnvkit.py batch "$BAM" \
    --reference "$REFERENCE" \
    --output-dir "$SAMPLE_DIR" \
    --scatter --diagram \
    --processes 4

echo "Finished CNVkit for $SAMPLE"

EOF

done

echo
echo "=== Submitted CNVkit jobs for all samples ==="
