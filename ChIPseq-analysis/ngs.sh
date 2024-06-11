#!/bin/bash
#$ -t 1-2

x=(
x
~/Sansom/Sansom_H3K27me3.txt
~/Sansom/Sansom_H3K4me3.txt
)

ngs.plot.r -G mm10 -R genebody -L 2000 -LEG 0 -C ${x[${SGE_TASK_ID}]} -O ${x[${SGE_TASK_ID}]}_1
ngs.plot.r -G mm10 -R tss -L 2000 -LEG 0 -C ${x[${SGE_TASK_ID}]} -O ${x[${SGE_TASK_ID}]}_2
