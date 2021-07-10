#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -N TE_count
#$ -o TE_count_output.txt
#$ -e TE_count_output.txt

## TE_HOME/data contains directories, each corresponding to a tissue. 
## Each tissue directory contains a group of aligned .bam files.

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/
cd $TE_HOME/data

## The first loop iterates through each tissue directory.

for TISSUE in *;do
  
  cd ${TISSUE}

  ## The second loop iterates through each .bam file for that tissue and runs
  ## TEcount, outputing a count table for gene and TE transcripts.
  
  for FILE in *.bam;do

    TEcount \
    -b ${FILE} \
    --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
    --TE $TE_HOME/index/annotations/TEtranscripts_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf \
    --project TEcount_${FILE}

  done

  ## After all the count tables for a tissue have been output, the following
  ## joins them into a single count table for each tissue replicate.
  
  declare -a file_array=()

  for TABLE in *cntTable; do
    file_array+=($TABLE)
  done

  counter=0
  while [ $((counter+1)) -lt ${#file_array[@]} ]; do

    if [ $counter -eq 0 ]; then

      join ${file_array[$counter]} ${file_array[$((counter+1))]} > output_$((counter+1)) 

    else

      join output_${counter} ${file_array[$((counter+1))]} > output_$((counter+1))

    fi 

    counter=$((counter+1))

  done

  cd ..
  
done

