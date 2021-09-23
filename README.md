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
    
## Endogenous retroviral elements in mTECs
ERE expression in medullary thymic epithelial cells.

## Dependencies
Python 3.9.6 
IPython 7.26.0
scipy 1.7.1
seaborn 0.11.1
matplotlib 3.4.2
matplotlib_venn 0.11.6
numpy 1.21.1
pandas 1.3.1
sklearn 0.24.2
statsmodels 0.12.2
