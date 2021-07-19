#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -N TE_count_individual
#$ -o TE_count_individual_output.txt
#$ -e TE_count_individual_output.txt

## TE_HOME/data contains directories, each corresponding to a tissue. 
## Each tissue directory contains a group of aligned .bam files. These
## files should be named in the following format:
##
## unique-identifier_tissue_batch_Aligned.out.bam
## e.g. SRR488685_ESC_UCSC_Aligned.out.bam
##
## TE_HOME/index/annotations contains two .gtf files:

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis

TISSUE=$1

echo $TISSUE

cd $TE_HOME/data/$TISSUE/

for FILE in *.bam;do

    TEcount \
    -b ${FILE} \
    --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
    --TE $TE_HOME/index/annotations/TEtranscripts_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf \
    --project TEcount_${FILE}

done

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
  
mv output_$((counter)) ../${TISSUE}.cntTable
rm output_*
  

