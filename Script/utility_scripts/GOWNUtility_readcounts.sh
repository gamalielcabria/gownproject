  GNU nano 2.9.8                                                                                 readcounts.bash                                                                                  Modified  

#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=readcount
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 
#SBATCH --time=48:00:00
####SBATCH --partition=bigmem
#SBATCH --mem=8G
############## Run your script ########


folderoffolders='/work/ebg_lab/eb/AIWIP_GOWN_amplicons'
for i in `ls -d ${folderoffolders}/reads_*`
do
    folder=${i##*/}
    output="readcounts_raw_${folder}.txt"

    for j in `ls ${i}/data/*fastq.gz`
    do
      	echo "${j} : `zcat ${j} | grep -c '^+'` "
    done > ${output}

done
