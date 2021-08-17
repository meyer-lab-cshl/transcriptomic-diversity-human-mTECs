#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=10G
#$ -N bam_to_bedgraph
#$ -o bam_to_bed_output.txt
#$ -e bam_to_bed_output.txt

cd $1

for FILE in *.bam; do
	
	NAME=`echo ${FILE} | cut -d'_' -f 1-3`
	
	samtools sort ${FILE} > ${NAME}.sorted.bam
	samtools index ${NAME}.sorted.bam ${NAME}.bam.bai
	
	bamCoverage \
	-b ${NAME}.sorted.bam \
	-of bedgraph \
	-o ${NAME}.bedgraph \
	--normalizeUsing CPM \
	--binSize 25 \
	--verbose

#	bedtools genomecov -ibam ${NAME}.sorted.bam -split -bg > ${NAME}.bedgraph
#	sort -k1,1 -k2,2n ${NAME}.bedgraph | bgzip > ${NAME}.bedgraph.gz
#	tabix -s 1 -b 2 -e 3 ${NAME}.bedgraph.gz

done
