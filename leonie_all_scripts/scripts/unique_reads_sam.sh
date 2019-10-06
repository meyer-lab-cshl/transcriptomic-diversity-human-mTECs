#!/bin/bash

dir=$filedir
name=$name

echo '1' 
/home/stroemic/software/miniconda2/envs/rnaseq/bin/samtools view -q 255 -h -b $dir/${name}_Aligned.out.sam > $dir/${name}_Aligned.out.unique.bam
echo '2'
/home/stroemic/software/miniconda2/envs/rnaseq/bin/samtools sort -o $dir/${name}_Aligned.out.unique.sorted.bam $dir/${name}_Aligned.out.unique.bam 
echo '3'
/home/stroemic/software/miniconda2/envs/rnaseq/bin/samtools index $dir/${name}_Aligned.out.unique.sorted.bam
