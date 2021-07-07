#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -N TE_count
#$ -o TE_count_output.txt
#$ -e TE_count_output.txt

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis
GTEX=/grid/meyer/home/jacarter/GTEx/Aligned_TE

cd $GTEX/

for TISSUE in *; do

  cd ${TISSUE}

  for FILE in GTEX*; do

    TEcount \
    -b ${FILE} \
    --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
    --TE $TE_HOME/index/annotations/TEtranscripts_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf \
    --project TEcount_~{TISSUE}_${FILE}

  done

done
