#!/bin/bash
#$ -t 1-7
#$ -l m_mem_free=8G
#$ -pe threads 4

#git clone https://github.com/BioinfoUNIBA/REDItools

#base=/grid/meyer/home/jacarter/TSS_RNAseq
#in=${base}/FASTQ
#out=${base}/FASTQ_merged
#cat ${in}/*221-hi*_1.fastq > ${out}/pt221_hi_1.fastq
#cat ${in}/*221-hi*_2.fastq > ${out}/pt221_hi_2.fastq
#cat ${in}/*221-lo*_1.fastq > ${out}/pt221_lo_1.fastq
#cat ${in}/*221-lo*_2.fastq > ${out}/pt221_lo_2.fastq
#cat ${in}/*226-hi*_1.fastq > ${out}/pt226_hi_1.fastq
#cat ${in}/*226-hi*_2.fastq > ${out}/pt226_hi_2.fastq
#cat ${in}/*226-lo*_1.fastq > ${out}/pt226_lo_1.fastq
#cat ${in}/*226-lo*_2.fastq > ${out}/pt226_lo_2.fastq
#cat ${in}/*Liver*_1.fastq > ${out}/Liver_1.fastq
#cat ${in}/*Liver*_2.fastq > ${out}/Liver_2.fastq
#cp ${in}/*Pt214-hi*_1.fastq ${out}/pt214_hi_1.fastq
#cp ${in}/*Pt214-hi*_2.fastq ${out}/pt214_hi_2.fastq
#cp ${in}/*Pt214-lo*_1.fastq ${out}/pt214_lo_1.fastq
#cp ${in}/*Pt214-lo*_2.fastq ${out}/pt214_lo_2.fastq

base=/grid/meyer/home/jacarter/TSS_RNAseq
in_dir=${base}/FASTQ_merged
files=($( ls ${in_dir}/*_1.fastq ))
subject="$(basename ${files[$SGE_TASK_ID-1]} _1.fastq)"

cd ${base}/FastQC
fastqc ${in_dir}/${subject}_1.fastq ${in_dir}/${subject}_2.fastq -o ${base}/FastQC

cd ${base}/FastP
fastp -i ${in_dir}/${subject}_1.fastq -I ${in_dir}/${subject}_2.fastq -o ${base}/FastP/${subject}_fastp_1.fastq -O ${base}/FastP/${subject}_fastp_2.fastq -q 25 -u 10 -l 50 -y -x -w 4

STAR_index=${HOME}/Reference_Genome/STAR_Index/hg38
aligned_dir=${base}/Aligned
STAR --runThreadN 4 \
    --runMode alignReads \
    --genomeDir ${STAR_index} \
    --genomeLoad NoSharedMemory \
    --outFileNamePrefix ${aligned_dir}/${subject}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --outFilterType BySJout \
    --outFilterMultimapNmax 1 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --readFilesIn ${base}/FastP/${subject}_fastp_1.fastq ${base}/FastP/${subject}_fastp_2.fastq
