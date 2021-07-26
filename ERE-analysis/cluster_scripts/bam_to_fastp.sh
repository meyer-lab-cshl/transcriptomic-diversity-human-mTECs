#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N bam_to_fastp
#$ -o bam_to_fastp_output.txt
#$ -e bam_to_fastp_output.txt

input=$1
new_name=`echo ${input} | cut -d'_' -f 1-3`

samtools sort -n -o $input $new_name.sorted.bam > $new_name.sorted.bam   

bedtools bamtofastq -i $new_name.sorted.bam \
                      -fq $new_name.R1.fq \
                      -fq2 $new_name.R2.fq

rm $new_name.sorted.bam

fastp -i $new_name.R1.fq -I $new_name.R2.fq \
-o $new_name.R1.fastp.fastq -O $new_name.R2.fastp.fastq \
-q 25 -u 10 -l 50 -y -x -w 4

rm $new_name.R1.fq
rm $new_name.R2.fq

mkdir fastp
mv *.fastp fastp
mv fastp.* fastp
mv bam_to_fastp_output.txt fastp
