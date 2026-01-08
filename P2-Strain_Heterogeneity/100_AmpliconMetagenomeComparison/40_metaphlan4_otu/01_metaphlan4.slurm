#!/bin/bash
# ---------------------------------------------------------------------
# An example PBS script for running a job on a compute cluster
# ---------------------------------------------------------------------
#SBATCH --job-name=metaphlan
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --time=72:00:00
##SBATCH --dependency=afterok:37497333
# ---------------------------------------------------------------------
echo "Starting run at: `date`"
# --------------------------------------------------------------------

# To run, install bbmap using `conda create --name bbmap -c bioconda bbmap
source ~/software/miniforge3/etc/profile.d/conda.sh 
conda activate metaphlan

## Setting up resources ##
input='/work/ebg_lab/eb/AIWIP_GOWN_metagenome/100_denitrificationindex/10_trimnfiltering/01_trimmed'
output='/work/ebg_lab/eb/AIWIP_GOWN_metagenome/200_MDSC507_metagenome_diversity/80_metaphlan_gc'
refDB='/home/gamaliellysandergl.c/database/metaphlan'
scratchdir=/scratch/${SLURM_JOBID}
tmp=/scratch/${SLURM_JOBID}/tmp
index='mpa_vJan25_CHOCOPhlAnSGB_202503'
env_dir='/home/gamaliellysandergl.c/software/miniforge3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/'

mkdir -p tmp
ls $tmp

## Running mmseqs2 ##
for i in `ls ${input}/*_1.clean.corrected.fastq.gz`;
do
filename=$(basename "$i" _1.clean.corrected.fastq.gz)
R1=${i}
R2=`echo $R1 | sed 's/_1/_2/g'`
Rs=`echo $R1 | sed 's/_1/_s/g'`
concat="/scatch/${SLURM_JOBID}/${filename}.concat.fastq.gz"

#cat $R1 $R2 $Rs > ${concat}

echo "Running ${filename}"

metaphlan --input_type fastq ${R1},${R2} \
-t rel_ab_w_read_stats \
--db_dir ${refDB} --index ${index} \
--mapout ${scratchdir}/${filename}.bowtie2.bz2 \
--tmp_dir ${tmp} --nproc ${SLURM_CPUS_PER_TASK} \
-o ${filename}_${index}.mp4.txt

# Modify taxa to match gtdb
${env_dir}/sgb_to_gtdb_profile.py -i ${filename}_${index}.mp4.txt -o ${filename}_${index}.mp4_gtdb.txt

done

# Merge mp4 output to one table 
absmerge='/home/gamaliellysandergl.c/database/metaphlan/merge_metaphlan_tables_abs.py'

echo "Merging profiles using ${absmerge}"

## Original SGB profiles
${env_dir}/merge_metaphlan_tables.py *.mp4.txt -o ${output}/abundance_table.mp4_merged_relab.tsv
${absmerge} *.mp4.txt -o ${output}/abundance_table.mp4_merged_abs.tsv

## GTDB profiles
${env_dir}/merge_metaphlan_tables.py --gtdb_profiles *.mp4_gtdb.txt -o ${output}/abundance_table.mp4_merged_gtdb_relab.tsv
${absmerge} --gtdb_profiles *.mp4_gtdb.txt -o ${output}/abundance_table.mp4_merged_gtdb_abs.tsv


# Calculate diversity Metrics
echo "Calculating Diversity"

Rscript ${env_dir}/calculate_diversity.R -f ${output}/abundance_table.mp4_merged_gtdb_abs.tsv -d beta -m bray-curtis -o ${output}/02_betadiversity
Rscript ${env_dir}/calculate_diversity.R -f ${output}/abundance_table.mp4_merged_gtdb_abs.tsv -d beta -m jaccard -o ${output}/02_betadiversity
Rscript ${env_dir}/calculate_diversity.R -f ${output}/abundance_table.mp4_merged_gtdb_abs.tsv -d beta -m aitchison -o ${output}/02_betadiversity

Rscript ${env_dir}/calculate_diversity.R -f ${output}/abundance_table.mp4_merged_gtdb_abs.tsv -d alpha -m richness -o ${output}/01_alphadiversity
Rscript ${env_dir}/calculate_diversity.R -f ${output}/abundance_table.mp4_merged_gtdb_abs.tsv -d alpha -m shannon -o ${output}/01_alphadiversity
Rscript ${env_dir}/calculate_diversity.R -f ${output}/abundance_table.mp4_merged_gtdb_abs.tsv -d alpha -m simpson -o ${output}/01_alphadiversity
Rscript ${env_dir}/calculate_diversity.R -f ${output}/abundance_table.mp4_merged_gtdb_abs.tsv -d alpha -m gini -o ${output}/01_alphadiversity

#richness, shannon, simpson, gini
