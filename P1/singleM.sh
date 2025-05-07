#!/bin/bash
## SCRIPT LOCATION = /work/ebg_lab/eb/strain_heterogeneity_project
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
readpath='/work/ebg_lab/eb/strain_heterogeneity_project/seqtk'

for i in `ls ${readpath}/*R1*fastq`
do
R2=`echo ${i} | sed 's/R1/R2/'`
R1=${i}
output_prefix="${R1%.*}.singlem"
output=${output_prefix%\.*}

#cat<<EOF
singlem pipe \
-1 ${R1} \
-2 ${R2} \
--taxonomic-profile-krona ${output}.krona.tsv \
--taxonomic-profile ${output}.taxa.txt \
--otu-table ${output}.otu.txt \
--threads 10

singlem summarise \
--input-taxonomic-profile ${output}.taxa.txt \
--output-taxonomic-profile-with-extras ${output}.summarise.wextras.tsv 

singlem summarise \
--input-taxonomic-profile ${output}.taxa.txt \
--output-species-by-site-relative-abundance-prefix ${output}_relabund

#singlem summarise \
#--input-otu-tables ${input}.otu.txt \
#--unifrac-by-otu ${input}_otu_unifrac-
#EOF

echo "SingleM run for ${metagenome} done at `date`"
done


#---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
