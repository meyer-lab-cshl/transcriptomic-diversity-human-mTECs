#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=30G
#$ -N TE_local_opt
#$ -o TE_local_optimization_output.txt
#$ -e TE_local_optimization_output.txt

cd $TE_HOME/data/STAR_optimization/mTEC-HI/                         

for FILE in *pt226*; do
  TElocal \
  -b ${FILE} \
  --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
  --TE $TE_HOME/index/annotations/TElocal_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf.locInd \
  --project TE_local_${FILE}

done
