#!/bin/bash
#$ -~TE_thymus/analysis/cluster/
#$ -pe threads 16
#$ -l m_mem_free=1G
#$ -N position_enrichment_permutation
#$ -o position_enrichment_permutation_output.txt
#$ -e position_enrichment_permutation_output.txt

Rscript position_enrichment_permutation.R
