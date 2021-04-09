#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=10G
#$ -N correlate_gene_and_TE_expression
#$ -o correlate_gene_and_TE_expression_output.txt
#$ -e correlate_gene_and_TE_expression_output.txt

module load EBModules
module load R/4.0.4-fosscuda-2020a
Rscript correlate_gene_and_TE_expression.R
