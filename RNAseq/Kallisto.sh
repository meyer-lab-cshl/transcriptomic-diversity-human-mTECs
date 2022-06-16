#!/bin/bash
#$ -t 1-7

base_dir=${HOME}/TSS_RNAseq
data_dir=${base_dir}/FastP
out_dir=${base_dir}/Kallisto/Quant
files=($( ls ${data_dir}/*_1.fastq ))
read1=${files[${SGE_TASK_ID}-1]}
read2=${read1/_1.fastq/_2.fastq}
filename="$(basename $read1 _1.fastq)"

mkdir -p ${out_dir}

index_dir=~/Reference_Genome/Kallisto_Index/
kallisto quant -i ${index_dir}transcripts.idx -o ${out_dir}/${filename} -b 100 ${read1} ${read2}
