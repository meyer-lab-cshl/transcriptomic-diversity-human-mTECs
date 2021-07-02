#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=3G
#$ -N index_genome
#$ -o index_genome_output.txt
#$ -e index_genome_output.txt

SAMPLE=['pt214_lo']

rule all:
  input: expand("{sample}...bam",
    sample=SAMPLE)

rule index_genome:
  input:
    gtf='index/human.GRCh38.gtf',
    fasta='index/human.GRCh38.genome.fa'
  output:
    'index/Genome'
  shell:
    """
      STAR \
      --sjdbOverhang 100 \
      --sjdbGTFfile  {input.gtf} \
      --runMode genomeGenerate \
      --genomeDir index \
      --genomeFastaFiles {input.fasta} \
      --runThreadN 4
    """

rule align:
  input:
    genome="index/Genome",
    fq1="{sample}_fastp_1.fastq",
    fq2="{sample}_fastp_2.fastq"
  output:
    "{sample}....bam"

