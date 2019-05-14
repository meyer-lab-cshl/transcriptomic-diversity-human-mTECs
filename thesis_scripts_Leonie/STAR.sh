#!/bin/bash

thread_number=$thread_number
genomedir=$genomedir
read1=$read1
#read2=$read2
outprefix=$outprefix

#/home/stroemic/software/miniconda2/envs/rnaseq/bin/STAR --runThreadN $thread_number --genomeDir $genomedir --readFilesIn $read1 $read2 --outFileNamePrefix $outprefix  
/home/stroemic/software/miniconda2/envs/rnaseq/bin/STAR --runThreadN $thread_number --genomeDir $genomedir --readFilesIn $read1 --outFileNamePrefix $outprefix  
