#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=10G
#$ -N count_overlaps
#$ -o count_overlaps_output.txt
#$ -e count_overlaps_output.txt

module load EBModules
module load R/4.0.4-fosscuda-2020a
Rscript count_overlaps.R




