#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=30G
#$ -N TE_small
#$ -o TE_small_output.txt
#$ -e TE_small_output.txt

cd $TE_HOME/data/

tesmall \
-f SRR7408176_parental_1.fastq \
-g hg38 \
-l TE_small_test
