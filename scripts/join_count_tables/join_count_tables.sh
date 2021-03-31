#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -N join_count_tables
#$ -o join_count_output.txt
#$ -e join_count_output.txt

cd /grid/meyer/home/mpeacey/TE_thymus/data/RNA_seq/testis_jason/count_tables

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

