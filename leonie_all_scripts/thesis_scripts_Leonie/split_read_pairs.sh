#!/bin/bash

xdir=$filedir
name=$name

echo '1' 
/home/stroemic/software/miniconda2/envs/rnaseq/bin/samtools view -h -f 0x40  $xdir/${name}_Aligned.out.unique.sorted.mbc.sam > $xdir/${name}_Aligned.out.unique.sorted.mbc._1.sam
echo '2'
/home/stroemic/software/miniconda2/envs/rnaseq/bin/samtools view -h -f 0x80  $xdir/${name}_Aligned.out.unique.sorted.mbc.sam > $xdir/${name}_Aligned.out.unique.sorted.mbc._2.sam
