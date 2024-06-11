#!/bin/bash
#$ -t 1-2

SRR=(
x
SRR1045005
SRR1045006
)

srr=${SRR[${SGE_TASK_ID}]}

picard MarkDuplicates \
-I ~/Sansom/Aligned/${srr}_Aligned.sortedByCoord.out.bam \
-O ~/Sansom/Rmdup/${srr}_Aligned.sortedByCoord.out.rmdup.bam \
-M ~/Sansom/Rmdup/Metrics/${srr}_dup_metrics.txt \
--REMOVE_DUPLICATES true
