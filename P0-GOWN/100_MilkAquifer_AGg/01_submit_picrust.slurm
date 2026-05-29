#!/bin/bash
#SBATCH --job-name=picrust2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=7
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu2025
#SBATCH --output=slurm-%x_%j.out
#SBATCH --error=slurm-%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gamaliellysandergl.c@ucalgary.ca

set -euo pipefail

echo "========================================"
echo "PICRUSt2 job started at: $(date)"
echo "========================================"

# Print full script for logging
sed 's/^/# /' "$0"

echo
echo "Runtime variables:"
echo "  JOB_ID:            ${SLURM_JOB_ID}"
echo "  JOB_NAME:          ${SLURM_JOB_NAME}"
echo "  NODES:             ${SLURM_JOB_NUM_NODES}"
echo "  CPUS_PER_TASK:     ${SLURM_CPUS_PER_TASK}"
echo "  MEMORY:            ${SLURM_MEM_PER_NODE:-$SLURM_MEM_PER_CPU}"
echo "  PARTITION:         ${SLURM_JOB_PARTITION}"
echo "  NODELIST:          ${SLURM_JOB_NODELIST}"
echo "  WORKDIR:           $(pwd)"
echo

# Conda init
source ~/software/miniforge3/etc/profile.d/conda.sh

# Activate PICRUSt2 environment
conda activate picrust2

# Inputs
FASTA="picrust_ASV.fasta"
OTU_TABLE="picrust_formatted_otu-table.tsv"
OUTDIR="picrust2_out"

echo "Input FASTA:     $FASTA"
echo "Input OTU table: $OTU_TABLE"
echo "Output folder:   $OUTDIR"
echo

echo "========================================"
echo "Running PICRUSt2 pipeline"
echo "========================================"

#picrust2_pipeline.py \
#  -s "$FASTA" \
#  -i "$OTU_TABLE" \
#  -o "$OUTDIR" \
#  -p "$SLURM_CPUS_PER_TASK"

echo
echo "========================================"
echo "Generating KO contribution table"
echo "========================================"

metagenome_pipeline.py \
  -i "$OTU_TABLE" \
  -f "$OUTDIR/combined_KO_predicted.tsv.gz" \
  -m "$OUTDIR/combined_marker_predicted_and_nsti.tsv.gz" \
  --strat_out \
  -o "$OUTDIR/KO_contrib_out"

echo
echo "========================================"
echo "Generating EC contribution table"
echo "========================================"

metagenome_pipeline.py \
  -i "$OTU_TABLE" \
  -f "$OUTDIR/combined_EC_predicted.tsv.gz" \
  -m "$OUTDIR/combined_marker_predicted_and_nsti.tsv.gz" \
  --strat_out \
  -o "$OUTDIR/EC_contrib_out"

echo
echo "========================================"
echo "Useful outputs"
echo "========================================"

echo "KO abundance table:"
echo "  $OUTDIR/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"

echo
echo "EC abundance table:"
echo "  $OUTDIR/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"

echo
echo "KO contribution table:"
echo "  $OUTDIR/KO_contrib_out/pred_metagenome_contrib.tsv.gz"

echo
echo "EC contribution table:"
echo "  $OUTDIR/EC_contrib_out/pred_metagenome_contrib.tsv.gz"

echo
echo "========================================"
echo "PICRUSt2 job finished at: $(date)"
echo "========================================"
