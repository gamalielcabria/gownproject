#!/bin/bash
#SBATCH --job-name=run260515_sweep ###combine_sweep
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200G
#SBATCH --time=5-00:00:00
#SBATCH --array=0-8
#SBATCH --partition=cpu2025
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --error=slurm-%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gamaliellysandergl.c@ucalgary.ca

set -euo pipefail

echo "========== SLURM JOB SCRIPT =========="
sed 's/^/# /' "$0"

echo "========== RUNTIME PARAMETERS =========="
echo "JOB_ID        : $SLURM_JOB_ID"
echo "ARRAY_JOB_ID  : ${SLURM_ARRAY_JOB_ID:-NA}"
echo "ARRAY_TASK_ID : ${SLURM_ARRAY_TASK_ID:-NA}"
echo "JOB_NAME      : $SLURM_JOB_NAME"
echo "NODES         : ${SLURM_JOB_NUM_NODES:-1}"
echo "CPUS_PER_TASK : $SLURM_CPUS_PER_TASK"
echo "MEMORY        : ${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-NA}}"
echo "PARTITION     : ${SLURM_JOB_PARTITION:-NA}"
echo "NODELIST      : ${SLURM_JOB_NODELIST:-NA}"
echo "WORKDIR       : $PWD"
echo "START TIME    : $(date)"
echo "======================================="

# --------------------------------------------------
# Paths
# --------------------------------------------------
RAW_BASE="/work/ebg_lab/eb/AIWIP_GOWN_amplicons/400_TimeSeries/02_renamedrawreads"
OUT_BASE="/work/ebg_lab/eb/AIWIP_GOWN_amplicons/400_TimeSeries/09_DADA2Sweep"
RSCRIPT="/work/ebg_lab/eb/AIWIP_GOWN_amplicons/400_TimeSeries/09_DADA2Sweep/01_DADA2_exploratory.R"

mkdir -p "$OUT_BASE"

echo "RAW_BASE : $RAW_BASE"
echo "OUT_BASE : $OUT_BASE"
echo "RSCRIPT  : $RSCRIPT"

# --------------------------------------------------
# Load conda
# --------------------------------------------------
source ~/software/miniforge3/etc/profile.d/conda.sh
conda activate dada2

echo "Using conda env:"
conda info --envs | grep '*'

# --------------------------------------------------
# Discover run* input folders
# --------------------------------------------------
mapfile -t RUN_DIRS < <(find "$RAW_BASE" -maxdepth 1 -mindepth 1 -type d -name 'run260515-GLC_24-25*' | sort)

if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: No run* directories found in $RAW_BASE" >&2
    exit 1
fi

echo "Found ${#RUN_DIRS[@]} run folders:"
printf '  %s\n' "${RUN_DIRS[@]}"

# --------------------------------------------------
# truncLen combinations
# --------------------------------------------------
TRUNC_F_VALS=(280 270 260)
TRUNC_R_VALS=(230 240 250)

N_RUNS=${#RUN_DIRS[@]}
N_COMBOS=$(( ${#TRUNC_F_VALS[@]} * ${#TRUNC_R_VALS[@]} ))   # 9
TOTAL_TASKS=$(( N_RUNS * N_COMBOS ))

TASK_ID=${SLURM_ARRAY_TASK_ID}

if (( TASK_ID < 0 || TASK_ID >= TOTAL_TASKS )); then
    echo "ERROR: TASK_ID $TASK_ID is out of range 0..$((TOTAL_TASKS-1))" >&2
    exit 1
fi

RUN_IDX=$(( TASK_ID / N_COMBOS ))
COMBO_IDX=$(( TASK_ID % N_COMBOS ))

F_IDX=$(( COMBO_IDX / ${#TRUNC_R_VALS[@]} ))
R_IDX=$(( COMBO_IDX % ${#TRUNC_R_VALS[@]} ))

INPUT_DIR="${RUN_DIRS[$RUN_IDX]}"
RUN_NAME="$(basename "$INPUT_DIR")"
TRUNC_F="${TRUNC_F_VALS[$F_IDX]}"
TRUNC_R="${TRUNC_R_VALS[$R_IDX]}"

OUTPUT_DIR="${OUT_BASE}/${RUN_NAME}/truncF${TRUNC_F}_truncR${TRUNC_R}"
mkdir -p "$OUTPUT_DIR"

echo "========== TASK MAPPING =========="
echo "TASK_ID    : $TASK_ID"
echo "RUN_IDX    : $RUN_IDX"
echo "COMBO_IDX  : $COMBO_IDX"
echo "RUN_NAME   : $RUN_NAME"
echo "INPUT_DIR  : $INPUT_DIR"
echo "TRUNC_F    : $TRUNC_F"
echo "TRUNC_R    : $TRUNC_R"
echo "OUTPUT_DIR : $OUTPUT_DIR"
echo "=================================="

Rscript "$RSCRIPT" "$INPUT_DIR" "$OUTPUT_DIR" "" "" "$TRUNC_F" "$TRUNC_R"

echo "Finished: $(date)"
