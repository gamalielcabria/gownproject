#!/bin/bash
## SCRIPT LOCATION = /work/ebg_lab/eb/strain_heterogeneity_project
####### Reserve computing resources #############
#SBATCH --job-name=SeqTK
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
####SBATCH --partition=bigmem
#SBATCH --mem=20G
############## Run your script ########


# Load SeqTK module if necessary which is in conda base env
source ~/software/init-conda
conda activate seqtk


# Define the folder where the raw reads are
#folder='/work/ebg_lab/eb/GOWN_samples_2019_2022/metagenomics_2022/qc'
folder='/work/ebg_lab/eb/strain_heterogeneity_project/seqtk'

# Define the file
filename=`ls ${folder}/*R1*.fastq`

# Get the prefix (left side of *.fastq)
prefix=`echo ${filename} | cut -f 1 -d 'R'`

# Subset 20M reads
i=1
for seed in 101 113 127;
do
seqtk sample -s ${seed} ${prefix}R1.fastq 20000000 > ${prefix}R1.20M.rep${i}.fastq
seqtk sample -s ${seed} ${prefix}R2.fastq 20000000 > ${prefix}R2.20M.rep${i}.fastq
((i++))
done

# Subset 10M reads
i=1
for seed in 137 149 163;
do
seqtk sample -s ${seed} ${prefix}R1.fastq 10000000 > ${prefix}R1.10M.rep${i}.fastq
seqtk sample -s ${seed} ${prefix}R2.fastq 10000000 > ${prefix}R2.10M.rep${i}.fastq
((i++))
done

# Subset 5M reads
i=1
for seed in 173 191 199; 
do
seqtk sample -s ${seed} ${prefix}R1.fastq 5000000 > ${prefix}R1.5M.rep${i}.fastq
seqtk sample -s ${seed} ${prefix}R2.fastq 5000000 > ${prefix}R2.5M.rep${i}.fastq
((i++))
done

# Subset 1M reads
i=1
for seed in 223 229 241;
do
seqtk sample -s ${seed} ${prefix}R1.fastq 1000000 > ${prefix}R1.1M.rep${i}.fastq
seqtk sample -s ${seed} ${prefix}R2.fastq 1000000 > ${prefix}R2.1M.rep${i}.fastq
((i++))
done

# Subset 500K reads
i=1
for seed in 251 263 277; 
do
seqtk sample ${prefix}R1.fastq 500000 > ${prefix}R1.500K.rep${i}.fastq
seqtk sample ${prefix}R2.fastq 500000 > ${prefix}R2.500K.rep${i}.fastq
((i++))
done

mkdir -p subset_fastq

mv *.fastq subset_fastq
