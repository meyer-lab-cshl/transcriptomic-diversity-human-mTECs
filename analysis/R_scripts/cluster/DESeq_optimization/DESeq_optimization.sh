#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=20G
#$ -N DESeq_optimization
#$ -o DESeq_optimization_output.txt
#$ -e DESeq_optimization_output.txt

module load EBModules
module load R/4.0.4-fosscuda-2020a
Rscript DESeq_optimization.R




