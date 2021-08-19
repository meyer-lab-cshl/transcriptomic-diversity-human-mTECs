#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -N RPKM_calculation
#$ -o RPKM_calculation_output.txt
#$ -e RPKM_calculation_output.txt

module load EBModules
module load R/4.0.4-fosscuda-2020a
Rscript RPKM_calculation.R
