#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=3G
#$ -N download_fastq_files
#$ -o download_fastq_files_output.txt
#$ -e download_fastq_files_output.txt

while read p; do
  fasterq-dump "$p" -p -e 15
done <acc.txt
