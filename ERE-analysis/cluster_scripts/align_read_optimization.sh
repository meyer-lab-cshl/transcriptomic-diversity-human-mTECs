#!/bin/bash
#$ -cwd
#$ -pe threads 4
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

outFilterMultimapNmax_array=(10 50 100 150 200)
winAnchorMultimapNmax_array=(20 100 200 300 400)

counter_1=0
for out_parameter in ${outFilterMultimapNmax_array[@]}; do
  
  counter_2=0
  for i in ${read_1_array[@]}; do

    STAR \
    --runThreadN 4 \
    --genomeDir $TE_HOME/index/ \
    --sjdbGTFfile $TE_HOME/index/human.GRCh38.gtf \
    --sjdbOverhang 100 \
    --readFilesIn $i ${read_2_array[$counter_2]}\
    --outSAMtype BAM Unsorted \
    --winAnchorMultimapNmax ${winAnchorMultimapNmax_array[$counter_1]} \
    --outFilterMultimapNmax $out_parameter \
    --outFileNamePrefix $TE_HOME/data/STAR_optimization/out_parameter_${out_parameter}/${out_parameter}_${i}_
    
    counter_2=$((counter_2+1))

  done

  counter_1=$((counter_1+1))
  
done

