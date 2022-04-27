# Transcriptome diversity in human medullary thymus epithelial cells

Using 5P-Cap-seq and RNAseq, we mapped transcription start regions and splicing
events in the thymus and compared them to corresponding maps for peripheral
tissues (FANTOM5 and GTEx data respectively).

This repository contains the data processing steps, from RNAseq reads to
transcription start site maps. It is divided by data into subdirectories.
Each subdirectory contains analysis organised in either a snakemake workflow
(with the workflow depicted as a DAG in the respective README.md) or
consecutively labelled scripts.

## 5P-Cap-seq-alignment
* Alignment and quality control of 5Pseq data from 5 human tissue donors

## 5P-Cap-seq-processing
*

## gene-lists
* gene lists derived from public datasets, and re-mapping to GRCh38
* estimating tissue specifc antigens based on tissue exclusivity index tau and
    22 tissues from GTEx
    
## Endogenous retroviral elements in mTECs
ERE expression in medullary thymic epithelial cells.

## Dependencies
- Python: python (v3.9.6), IPython (v7.26.0), scipy (v1.7.1), seaborn (v0.11.1), matplotlib (v3.4.2), matplotlib_venn (v0.11.6), numpy (v1.21.1), pandas (v1.3.1),  sklearn (v0.24.2), statsmodels (v0.12.2)
- R: R (v4.0.3), CAGEr (v1.32), limma (v3.46), Sleuth (v0.30.0), biomaRt (v2.46.3), rMATS (v4.1.1)
- sequence analysis: umi_tools (v1.1), fastq_screen (v0.14.0), STAR (v2.7.2b), samtools (v1.11), picard (v2.18.20), multiqc (v1.9), paraclu (v9), bedops (v2.4.38), bedtools (v2.29.2), HOMER (v4.11.1), sambamba (v0.8), fastp (v0.11.8), Kallisto (v0.4.6)
