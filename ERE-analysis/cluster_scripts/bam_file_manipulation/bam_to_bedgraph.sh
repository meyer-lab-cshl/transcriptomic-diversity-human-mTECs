#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=16G
#$ -N bam_to_bedgraph
#$ -o bam_to_bed_output.txt
#$ -e bam_to_bed_output.txt

cd $1
mkdir bedgraph

for FILE in *Aligned.out.bam; do
	
	NAME=`echo ${FILE} | cut -d'_' -f 1-3`
	
	#samtools sort ${FILE} > ${NAME}.sorted.bam
	#samtools index ${NAME}.sorted.bam ${NAME}.sorted.bam.bai
	
	bamCoverage \
	-b ${NAME}.sorted.bam \
	-of bedgraph \
	-o bedgraph/${NAME}.bedgraph \
	--normalizeUsing CPM \
	--binSize 25 \
	--verbose \
	-p 4

#	bedtools genomecov -ibam ${NAME}.sorted.bam -split -bg > ${NAME}.bedgraph
#	sort -k1,1 -k2,2n ${NAME}.bedgraph | bgzip > ${NAME}.bedgraph.gz
#	tabix -s 1 -b 2 -e 3 ${NAME}.bedgraph.gz

done
