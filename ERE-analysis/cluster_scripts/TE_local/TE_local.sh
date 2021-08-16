#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -N TE_local
#$ -o TE_local_output.txt
#$ -e TE_local_output.txt

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis

cd $TE_HOME/data/${1}

TElocal \
-b ${2} \
--GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
--TE $TE_HOME/index/annotations/TElocal_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf.locInd \
--project TE_local_${2}




