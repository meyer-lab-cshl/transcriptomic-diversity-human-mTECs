#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -N rename_fastp
#$ -o rename_fastp_output.txt
#$ -e rename_fastp_output.txt

cd $1

for FILE in *.fastp;do

	mv ${FILE} ${FILE}.fastq

done

