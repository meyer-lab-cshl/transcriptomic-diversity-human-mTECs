#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=4G
#$ -N TE_transcripts
#$ -o TE_transcripts_output.txt
#$ -e TE_transcripts_output.txt

cd $TE_HOME/data

TEtranscripts \
--format BAM \
-t pt214_hi_fastp_1.fastq_Aligned.out.bam pt221_hi_fastp_1.fastq_Aligned.out.bam pt226_hi_fastp_1.fastq_Aligned.out.bam \
-c pt214_lo_fastp_1.fastq_Aligned.out.bam pt221_lo_fastp_1.fastq_Aligned.out.bam pt226_lo_fastp_1.fastq_Aligned.out.bam \
--GTF $TE_HOME/index/human.GRCh38.gtf \
--TE $TE_HOME/index/GRCh38_Ensembl_rmsk_TE.gtf \
--mode multi \
--project hi_vs_lo

