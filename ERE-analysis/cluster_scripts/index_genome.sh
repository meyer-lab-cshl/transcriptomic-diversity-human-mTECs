#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=3G
#$ -N index_genome
#$ -o index_genome_output.txt
#$ -e index_genome_output.txt

STAR \
--sjdbOverhang 100 \
--sjdbGTFfile $TE_HOME/index/human.GRCh38.gtf \
--runMode genomeGenerate \
--genomeDir $TE_HOME/index \
--genomeFastaFiles $TE_HOME/index/human.GRCh38.genome.fa \
--runThreadN 4
