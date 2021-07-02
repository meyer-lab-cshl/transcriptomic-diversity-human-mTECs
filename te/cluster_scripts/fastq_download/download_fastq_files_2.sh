#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=5G
#$ -N download_fastq_files
#$ -o download_fastq_files_2_output.txt
#$ -e download_fastq_files_2_output.txt

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/SRR488684/SRR488684_1.fastq.gz -o SRR488684_polyA_RNA_sequencing_of_H1_cells_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/SRR488684/SRR488684_2.fastq.gz -o SRR488684_polyA_RNA_sequencing_of_H1_cells_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/SRR488685/SRR488685_1.fastq.gz -o SRR488685_polyA_RNA_sequencing_of_H1_cells_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR488/SRR488685/SRR488685_2.fastq.gz -o SRR488685_polyA_RNA_sequencing_of_H1_cells_2.fastq.gz
