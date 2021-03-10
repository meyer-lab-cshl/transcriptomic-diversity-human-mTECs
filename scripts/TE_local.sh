#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=30G
#$ -N TE_local
#$ -o TE_local_output.txt
#$ -e TE_local_output.txt

cd $TE_HOME/data

TEcount \
-b pt214_hi_fastp_1.fastq_Aligned.out.bam \
--GTF $TE_HOME/index/human.GRCh38.gtf \
--TE $TE_HOME/index/GRCm38_Ensembl_rmsk_TE.gtf.locInd \
--project TE_local_test
