#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=30G
#$ -N TE_local
#$ -o TE_local_output.txt
#$ -e TE_local_output.txt

cd $TE_HOME/data/STAR_optimization/

for FILE in *_Aligned.out.bam; do
  TElocal \
  -b ${FILE} \
  --GTF $TE_HOME/index/human.GRCh38.gtf \
  --TE $TE_HOME/index/hg38_rmsk_TE.gtf.locInd \
  --project TE_local_${FILE}

done
