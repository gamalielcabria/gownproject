#!/bin/bash
# ---------------------------------------------------------------------
# An example PBS script for running a job on a compute cluster
# ---------------------------------------------------------------------
#SBATCH --job-name=singlem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=150gb
#SBATCH --time=7-00:00:00

# ---------------------------------------------------------------------
echo "Starting run at: `date`"
# ---------------------------------------------------------------------

source ~/software/miniforge3/etc/profile.d/conda.sh 
conda activate singlem

export SINGLEM_METAPACKAGE_PATH='/home/gamaliellysandergl.c/software/singlem/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'
#readpath='/home/gamaliellysandergl.c/p2_workstation/reads'

for i in `ls *taxa.txt`
do
#metagenome=${i##/*}
#R1=`ls ${readpath}/${metagenome}/*.qc.R1.fastq`
#R2=`ls ${readpath}/${metagenome}/*.qc.R2.fastq`
#output="${metagenome}_singlem"

#cat<<EOF
#singlem pipe \
#-1 ${R1} \
#-2 ${R2} \
#--taxonomic-profile-krona ${output}.krona.tsv \
#--taxonomic-profile ${output}.taxa.txt \
#--otu-table ${output}.otu.txt \
#--threads 10

input=${i%%\.*}

#cat<<EOF
singlem summarise \
--input-taxonomic-profile ${input}.taxa.txt \
--output-taxonomic-profile-with-extras ${input}.summarise.wextras.tsv 

singlem summarise \
--input-taxonomic-profile ${input}.taxa.txt \
--output-species-by-site-relative-abundance-prefix ${input}_relabund

#singlem summarise \
#--input-otu-tables ${input}.otu.txt \
#--unifrac-by-otu ${input}_otu_unifrac-

#EOF



echo "SingleM run for ${metagenome} done at `date`"
done


#---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
