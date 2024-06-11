#!/bin/bash
#$ -l m_mem_free=5G
#$ -pe threads 4

base=${HOME}/TSS_RNAseq

gtf=${HOME}/Reference_Genome/Genome_and_Annotation/hg38/hg38.gtf
star_index=${HOME}/Reference_Genome/STAR_Index/hg38
s1=${base}/Code/S1.txt
s2=${base}/Code/S2.txt
out=${base}/RMATS
temp=/grid/meyer/home/jacarter/TSS_RNAseq/RMATS/temp

${HOME}/ToolKit/rmats_turbo_v4_1_1/run_rmats --b1 ${s1} --b2 ${s2} --gtf ${gtf} --bi ${star_index} -t paired --readLength 51 --variable-read-length --nthread 4 --od ${out} --tmp ${temp}
