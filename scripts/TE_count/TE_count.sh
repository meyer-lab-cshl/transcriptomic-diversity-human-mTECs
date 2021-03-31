#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G
#$ -N TE_count_T
#$ -o TE_count_output.txt
#$ -e TE_count_output.txt

tissue=testis_jason

cd $TE_HOME/data/RNA_seq/${tissue}/

for FILE in GTEX*; do

  TEcount \
  -b ${FILE} \
  --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
  --TE $TE_HOME/index/annotations/TEtranscripts_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf \
  --project TE_count_${FILE}

done
