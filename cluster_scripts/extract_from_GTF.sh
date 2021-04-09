#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=3G
#$ -N extract_from_GTF
#$ -o extract_from_GTF_output.txt
#$ -e extract_from_GTF_output.txt

cd $TE_HOME/index

cat human.GRCh38.gtf | \
awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' | \
sed 's/gene_id "//' | \
sed 's/gene_id "//' | \
sed 's/gene_type "//'| \
sed 's/gene_name "//' | \
sed 's/"//g' | \
awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | \
sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" \
> gencode.v38_gene_annotation_table.txt
