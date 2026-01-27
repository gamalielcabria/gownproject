#!/bin/bash
#SBATCH --job-name=qiime2_dada2_sweep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=qiime2_sweep_%j.out
#SBATCH --error=qiime2_sweep_%j.err
#SBATCH --partition=cpu2023

set -euo pipefail

# -------------------------
# Inputs
# -------------------------
MANIFEST="manifest.tsv"
METADATA="sample-metadata.tsv"
CLASSIFIER="gtdb_classifier_r220.qza"
MANIFEST_FORMAT="PairedEndFastqManifestPhred33V2"

# Base output directory for the sweep
BASE_OUTDIR="21_qiime2_asv_gtdb_sweep_24-25"
mkdir -p "$BASE_OUTDIR"

# -------------------------
# DADA2 fixed params (same across all runs)
# -------------------------
TRIM_LEFT_F=0
TRIM_LEFT_R=0

MAX_EE_F=2
MAX_EE_R=2
TRUNC_Q=2
CHIMERA_METHOD="consensus"   # consensus|pooled|none
MIN_OVERLAP=12
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# -------------------------
# Define your trunc-len tests
# -------------------------
# Option A: explicit pairs (recommended)
# Put the combos you actually want to compare.
PAIRS=(
  "280 260"
  "270 250"
  "260 240"
  "250 230"
)

# Option B: grid search (uncomment if you want a grid)
# F_LIST=(300 290 280 270 260 250)
# R_LIST=(280 270 260 250 240 230)

# -------------------------
# Environment
# -------------------------
source /home/gamaliellysandergl.c/software/init-conda
conda activate qiime2-amplicon-2024.5

export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

# -------------------------
# Checks
# -------------------------
for f in "$MANIFEST" "$METADATA" "$CLASSIFIER"; do
  [[ -f "$f" ]] || { echo "ERROR: missing file: $f" >&2; exit 1; }
done

echo "QIIME2 version:"
qiime --version

# -------------------------
# 1) Import reads ONCE
# -------------------------
DEMUX_QZA="$BASE_OUTDIR/demux.qza"
DEMUX_QZV="$BASE_OUTDIR/demux.qzv"

if [[ ! -f "$DEMUX_QZA" ]]; then
  echo "Importing reads -> $DEMUX_QZA"
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --input-format "$MANIFEST_FORMAT" \
    --output-path "$DEMUX_QZA"
else
  echo "Found existing $DEMUX_QZA (skipping import)"
fi

if [[ ! -f "$DEMUX_QZV" ]]; then
  qiime demux summarize \
    --i-data "$DEMUX_QZA" \
    --o-visualization "$DEMUX_QZV"
fi

# -------------------------
# Summary file header
# -------------------------
SUMMARY="$BASE_OUTDIR/trunc_sweep_summary.tsv"
echo -e "trunc_f\ttrunc_r\tinput\tfiltered\tdenoised\tmerged\tnonchim\tasvs" > "$SUMMARY"

# Helper: export denoising-stats to tsv and sum across samples
sum_denoise_stats() {
  local denoise_qza="$1"
  local tmpdir="$2"

  mkdir -p "$tmpdir"
  qiime tools export --input-path "$denoise_qza" --output-path "$tmpdir" >/dev/null

  # exported file is usually: stats.tsv
  local tsv="$tmpdir/stats.tsv"
  [[ -f "$tsv" ]] || { echo "ERROR: cannot find $tsv after exporting denoising stats" >&2; exit 1; }

  # Sum numeric columns across all samples (skip header)
  # Columns in stats.tsv typically:
  # sample-id, input, filtered, denoised, merged, non-chimeric
  awk -F'\t' 'NR>1{
    in+=$2; filt+=$3; den+=$4; mer+=$5; non+=$6
  } END{
    printf "%d\t%d\t%d\t%d\t%d", in, filt, den, mer, non
  }' "$tsv"
}

# Helper: count ASVs (features) from rep-seqs export FASTA
count_asvs_from_rep_fasta() {
  local rep_qza="$1"
  local tmpdir="$2"

  mkdir -p "$tmpdir"
  qiime tools export --input-path "$rep_qza" --output-path "$tmpdir" >/dev/null

  local fasta="$tmpdir/dna-sequences.fasta"
  [[ -f "$fasta" ]] || { echo "ERROR: cannot find $fasta after exporting rep seqs" >&2; exit 1; }

  grep -c '^>' "$fasta"
}

# -------------------------
# 2) Run sweep
# -------------------------
run_one() {
  local TRUNC_LEN_F="$1"
  local TRUNC_LEN_R="$2"

  local TAG="truncF${TRUNC_LEN_F}_truncR${TRUNC_LEN_R}"
  local OUTDIR="$BASE_OUTDIR/$TAG"
  mkdir -p "$OUTDIR"

  echo
  echo "=============================="
  echo "Running DADA2: $TAG"
  echo "Output: $OUTDIR"
  echo "=============================="

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$DEMUX_QZA" \
    --p-trim-left-f "$TRIM_LEFT_F" \
    --p-trim-left-r "$TRIM_LEFT_R" \
    --p-trunc-len-f "$TRUNC_LEN_F" \
    --p-trunc-len-r "$TRUNC_LEN_R" \
    --p-max-ee-f "$MAX_EE_F" \
    --p-max-ee-r "$MAX_EE_R" \
    --p-trunc-q "$TRUNC_Q" \
    --p-chimera-method "$CHIMERA_METHOD" \
    --p-min-overlap "$MIN_OVERLAP" \
    --p-n-threads 0 \
    --o-table "$OUTDIR/asv-table.qza" \
    --o-representative-sequences "$OUTDIR/asv-rep-seqs.qza" \
    --o-denoising-stats "$OUTDIR/denoising-stats.qza"

  # Export + summarize stats
  local sums
  sums=$(sum_denoise_stats "$OUTDIR/denoising-stats.qza" "$OUTDIR/_export_denoise_stats")

  # Count ASVs
  local asvs
  asvs=$(count_asvs_from_rep_fasta "$OUTDIR/asv-rep-seqs.qza" "$OUTDIR/_export_rep_seqs")

  echo -e "${TRUNC_LEN_F}\t${TRUNC_LEN_R}\t${sums}\t${asvs}" >> "$SUMMARY"
}

# Run explicit pairs:
for pair in "${PAIRS[@]}"; do
  run_one $pair
done

# Or run grid search (uncomment if using F_LIST/R_LIST):
# for f in "${F_LIST[@]}"; do
#   for r in "${R_LIST[@]}"; do
#     run_one "$f" "$r"
#   done
# done

echo
echo "DONE."
echo "Summary table:"
echo "  $SUMMARY"
echo "Per-run outputs under:"
echo "  $BASE_OUTDIR/truncF*_truncR*/"
