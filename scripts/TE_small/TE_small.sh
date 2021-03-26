#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=30G
#$ -N TE_small
#$ -o TE_small_output.txt
#$ -e TE_small_output.txt

cd $TE_HOME/data/5Pseq/raw_fastq/

tesmall \
-f pt212-hi_1.paired.fastq \
-g hg38 \
-l test
