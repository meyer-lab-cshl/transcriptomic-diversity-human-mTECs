#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=32G
#$ -N align_reads
#$ -o align_reads_output.txt
#$ -e align_reads_output.txt

tissue=ESCs

cd $TE_HOME/data/RNA_seq/${tissue}/raw_fastq/

declare -a read_1_array=()
declare -a read_2_array=()

for READ_1 in *_1.fastq; do
  read_1_array+=($READ_1)
done

for READ_2 in *_2.fastq; do
  read_2_array+=($READ_2)
done

## Iterate through each array and run fastp/STAR on each pair of reads

counter=0
for i in ${read_1_array[@]}; do

  fastp \
  -i $TE_HOME/data/RNA_seq/${tissue}/raw_fastq/${i} \
  -I $TE_HOME/data/RNA_seq/${tissue}/raw_fastq/${read_2_array[$counter]} \
  -o $TE_HOME/data/RNA_seq/${tissue}/fastp/${i}.fastp \
  -O $TE_HOME/data/RNA_seq/${tissue}/fastp/${read_2_array[$counter]}.fastp \
  -q 25 -u 10 -l 50 -y -x -w 4
  
  STAR \
  --runThreadN 4 \
  --genomeDir $TE_HOME/index/ \
  --sjdbGTFfile $TE_HOME/index/annotations/human.GRCh38.gtf \
  --sjdbOverhang 100 \
  --readFilesIn $TE_HOME/data/RNA_seq/${tissue}/fastp/${i}.fastp \
  $TE_HOME/data/RNA_seq/${tissue}/fastp/${read_2_array[$counter]}.fastp \
  --outSAMtype BAM Unsorted \
  --winAnchorMultimapNmax 200 \
  --outFilterMultimapNmax 100 \
  --outFileNamePrefix $TE_HOME/data/RNA_seq/${tissue}/bam_files/${i}_
  
  counter=$((counter+1))

done

