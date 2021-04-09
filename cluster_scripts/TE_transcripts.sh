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
-t $TE_HOME/data/GTEX_test/Not_Sun_Exposed/GTEX-15ER7-0626-SM-6PANF_Aligned.out.bam $TE_HOME/data/GTEX_test/Not_Sun_Exposed/GTEX-1AX9J-1326-SM-731BJ_Aligned.out.bam \
-c $TE_HOME/data/GTEX_test/Brain/GTEX-132Q8-0011-R7b-SM-5N9F1_Aligned.out.bam $TE_HOME/data/GTEX_test/Brain/GTEX-13X6K-0011-R7b-SM-5P9K7_Aligned.out.bam \
--GTF $TE_HOME/index/human.GRCh38.gtf \
--TE $TE_HOME/index/GRCh38_GENCODE_rmsk_TE.gtf \
--mode multi \
--project GTEX_test
