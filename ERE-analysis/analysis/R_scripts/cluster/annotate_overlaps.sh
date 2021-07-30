#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -N annotate_overlaps
#$ -o annotate_overlaps_output.txt
#$ -e annotate_overlaps__output.txt

module load EBModules
module load R/4.0.4-fosscuda-2020a
Rscript annotate_overlaps.R
