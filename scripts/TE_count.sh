#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=8G
#$ -N TE_count
#$ -o TE_count_output.txt
#$ -e TE_count_output.txt

cd $TE_HOME/data

TEcount \
-b pt214_hi_fastp_1.fastq_Aligned.out.bam \
--GTF $TE_HOME/index/human.GRCh38.gtf \
--TE $TE_HOME/index/GRCh38_GENCODE_rmsk_TE.gtf \
--project TE_count_test

