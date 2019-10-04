#!/bin/bash

xdir=$filedir
name=$name
strand=$strand

echo '1' 
/home/stroemic/software/miniconda2/envs/rnaseq/bin/samtools view -h -b $xdir/${name}_Aligned.out.unique.sorted.mbc.collapsed_1.${strand}.sam > $xdir/${name}_Aligned.out.unique.sorted.mbc.collapsed_1.${strand}.bam
echo '2'
/home/stroemic/software/miniconda2/envs/rnaseq/bin/samtools sort -o $xdir/${name}_Aligned.out.unique.sorted.mbc.collapsed_1.${strand}.sorted.bam $xdir/${name}_Aligned.out.unique.sorted.mbc.collapsed_1.${strand}.bam 
echo '3'
/home/stroemic/software/miniconda2/envs/rnaseq/bin/samtools index $xdir/${name}_Aligned.out.unique.sorted.mbc.collapsed_1.${strand}.sorted.bam
