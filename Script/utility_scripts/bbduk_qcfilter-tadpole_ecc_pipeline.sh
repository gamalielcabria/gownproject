#!/bin/bash
# BBDUK Cleanup in Quadrusaurus

resources='/home/charlie/miniforge3/envs/bbmap_env/opt/bbmap-39.01-1/resources'

### 02 - Adapter Removal ###
for i in `ls *1.fq.gz`; 
do 

n=${i%_*}


bbduk.sh t=7 in1=${n}_1.fq.gz in2=${n}_2.fq.gz \
out1=${n}_1.trim.fq.gz out2=${n}_2.trim.fq.gz outs=${n}_s.trim.fq.gz \
ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 \
ref=${resources}/adapters.fa ordered ow=t stats=${n}.adapter.txt

echo "Adapter removed: ${i} \n"

done

### 03 - Remove contaminants  ###
for i in `ls *1.trim.fq.gz`; 
do 

n=${i%_*}


bbduk.sh t=7 in1=${n}_1.trim.fq.gz in2=${n}_2.trim.fq.gz \
out1=${n}_1.decon.fq.gz out2=${n}_2.decon.fq.gz outs=${n}_s1.decon.fq.gz \
k=31 ref=${resources}/sequencing_artifacts.fa.gz,${resources}/phix174_ill.ref.fa.gz ordered cardinality ow=t stats=${n}.decon.txt


bbduk.sh t=7 in=${n}_s.trim.fq.gz out1=${n}_s2.decon.fq.gz \
k=31 ref=${resources}/sequencing_artifacts.fa.gz,${resources}/phix174_ill.ref.fa.gz ordered cardinality ow=t stats=${n}.s.decon.txt

cat ${n}_s1.decon.fq.gz ${n}_s2.decon.fq.gz > ${n}_s.decon.fq.gz
rm ${n}_s1.decon.fq.gz ${n}_s2.decon.fq.gz

echo "Contaminants removed: ${i} \n"

done

### 04 - Remove Low Quality Reads  ###
for i in `ls *1.decon.fq.gz`; 
do 

n=${i%_*}


bbduk.sh t=7 in1=${n}_1.decon.fq.gz in2=${n}_2.decon.fq.gz \
out1=${n}_1.clean.fq.gz out2=${n}_2.clean.fq.gz outs=${n}_s1.clean.fq.gz \
trimq=20 maq=10 minlen=70 ordered maxns=0 ow=t


bbduk.sh t=7 in=${n}_s.decon.fq.gz out1=${n}_s2.clean.fq.gz \
trimq=20 maq=10 minlen=70 ordered maxns=0 ow=t

cat ${n}_s1.clean.fq.gz ${n}_s2.clean.fq.gz > ${n}_s.clean.fq.gz
rm ${n}_s1.clean.fq.gz ${n}_s2.clean.fq.gz

echo "Contaminants removed: ${i} \n"

done

### 5 - Error Correction  ###
for i in `ls *1.clean.fq.gz`; 
do 

n=${i%_*}

tadpole.sh t=7 in1=${n}_1.clean.fq.gz in2=${n}_2.clean.fq.gz out1=${n}_1.clean.correct.fq.gz out2=${n}_2.clean.correct.fq.gz mode=correct k=62 -Xmx100g prealloc=t |& tee -a ${n}_read_correct.out.txt


done

