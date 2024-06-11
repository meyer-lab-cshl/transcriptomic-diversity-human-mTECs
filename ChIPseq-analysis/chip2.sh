#!/bin/bash
#$ -t 1-6

SRR=(
x
SRR1045003
SRR1045004
SRR1045007
SRR1045008
)

srr=${SRR[${SGE_TASK_ID-1}]}

data_dir=/grid/meyer/home/jacarter/Sansom
out_dir=${data_dir}/Aligned

cd ${data_dir}/Raw
fastq-dump --split-files $srr

read1=${data_dir}/Raw/${srr}_1.fastq
read2=${data_dir}/Raw/${srr}_2.fastq

read1_fastp=${data_dir}/FASTP/${srr}_1.paired.fastp.fastq
read2_fastp=${data_dir}/FASTP/${srr}_2.paired.fastp.fastq

fastp -i ${read1} -I ${read2} -o ${read1_fastp} -O ${read2_fastp}

STAR_index=${HOME}/Reference_Genome/STAR_Index/mm10

STAR --runThreadN 8 \
    --runMode alignReads \
    --genomeDir ${STAR_index} \
    --readFilesIn ${read1_fastp} ${read2_fastp} \
    --outReadsUnmapped Fastq \
    --outSAMattrRGline ID:${srr} SM:${srr} \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix ${out_dir}/${srr}_
