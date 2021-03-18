#!/bin/bash
#$ -cwd
#$ -pe threads 16
#$ -l m_mem_free=1G
#$ -N position_enrichment_permutation
#$ -o position_enrichment_permutation_output.txt
#$ -e position_enrichment_permutation_output.txt

module load EBModules
module load R/4.0.4-fosscuda-2020a
Rscript position_enrichment_permutation.R
