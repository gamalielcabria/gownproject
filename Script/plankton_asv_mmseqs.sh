#!/bin/bash
# ---------------------------------------------------------------------
# An example PBS script for running a job on a compute cluster
# ---------------------------------------------------------------------
#SBATCH --job-name=mmseq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=24G
#SBATCH --time=72:00:00
# ---------------------------------------------------------------------
echo "Starting run at: `date`"
# --------------------------------------------------------------------
# To run, install bbmap using `conda create --name bbmap -c bioconda bbmap
source ~/software/miniforge3/etc/profile.d/conda.sh 
conda activate mmseqs2

query_path="OldConsensusSeqsForTest2025Dec08.fasta"
target_path="SOZoops18SOperon-ReferenceSeqs.fasta"
#taxon_path="insert_path_to/taxonomy_results.csv"
query="$(basename "${query_path%.*}")"
target="$(basename "${target_path%.*}")"

SA
mmseqs easy-search $query $target ${query}_v_${target}.output.m8 tmp --cov-mode 2 -c 0.97 --min-seq-id 0.98 --search-type 3

#Rscript
#Rscript plankton_asv_match.R $taxon_path ${query}_v_${target}.output.m8

# - --------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------