#!/bin/bash
#$ -cwd
#$ -pe threads 12
#$ -l m_mem_free=32G
#$ -N align_reads
#$ -o align_reads_output.txt
#$ -e align_reads_output.txt

cd $TE_HOME/data/FastP

## Generate arrays to hold paired reads

declare -a read_1_array=()
declare -a read_2_array=()

for READ_1 in *_1.fastq; do
  read_1_array+=($READ_1)
done

for READ_2 in *_2.fastq; do
  read_2_array+=($READ_2)
done

## Iterate through each array and run STAR on each pair of reads

counter=0
for i in ${read_1_array[@]}; do

  STAR \
  --runThreadN 12 \
  --genomeDir $TE_HOME/index/ \
  --sjdbGTFfile $TE_HOME/index/human.GRCh38.gtf \
  --sjdbOverhang 100 \
  --readFilesIn $i ${read_2_array[$counter]}\
  --outSAMtype BAM Unsorted \
  --winAnchorMultimapNmax 200 \
  --outFilterMultimapNmax 100 \
  --outFileNamePrefix $TE_HOME/data/${i}_
  
  counter=$((counter+1))

done
