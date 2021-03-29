#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -N TE_count
#$ -o TE_count_output.txt
#$ -e TE_count_output.txt

tissue=testis_jason

cd $TE_HOME/data/RNA_seq/${tissue}/test

counter=0
for FILE in GTEX*; do

  TEcount \
  -b ${FILE} \
  --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
  --TE $TE_HOME/index/annotations/TEtranscripts_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf \
  --project TE_count_${FILE}

  if [${counter} -eq 0]
  then
    TE_count_${FILE} > count_table
  else
    join count_table TE_count_${FILE} > count_table
  fi

  counter=$((counter+1))

done
