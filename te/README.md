# Increased mis-initiation and non-promiscuous splicing in human medullary thymus epithelial cells

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

=======
## Endogenous retroviral elements in mTECs
Transposable element expression in medullary thymic epithelial cells.

## Workflow:

1. **Pre-processing with FastP**

2. **Align reads with STAR**
    + GENCODE GRCh38 gene annotation
    + sjdbOverhang: 100
    + winAnchorMultimapNMax: 200
    + outFilterMultimapNmax: 100
    
3. **Quantify gene and TE expression with TEcounts or TElocal**
    + GENCODE GRCh38 Repeatmasker TE annotation (downloaded from http://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/prebuilt_indices/)
    + GENCODE GRCh38 gene annotation
    
4. **Differential expression analysis with DESEq2**
    + Pre-filtering: minimum of 2 normalized reads per gene/TE
    + FDR: 0.1
    + Independent filtering: off
