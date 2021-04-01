#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -N TE_count_ESCs
#$ -o TE_count_2_output.txt
#$ -e TE_count_2_output.txt

tissue=ESCs 

cd $TE_HOME/data/RNA_seq/${tissue}/bam_files 

for FILE in *bam; do

  TEcount \
  -b ${FILE} \
  --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
  --TE $TE_HOME/index/annotations/TEtranscripts_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf \
  --project TE_count_${FILE}

done

declare -a file_array=()

for FILE in *cntTable; do
  file_array+=($FILE)
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


