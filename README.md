# Transcriptomic diversity in human medullary thymus epithelial cells

## Study
The induction of central T-cell tolerance in the thymus depends on the  presentation of peripheral self-epitopes by medullary thymic epithelial cells (mTECs), enabled by a process known as promiscuous gene expression (pGE). Transcriptome diversity generated during pGE has many contributors, including non-canonical transcription initiation, alternative splicing and expression of endogenous retroelements (EREs). However, their significance and regulation are poorly understood in the healthy human thymus. Here, we mapped the expression of genome-wide transcripts in immature and mature human mTECs using high-throughput 5'Cap and RNA sequencing. Overall, 96\% of protein coding genes were represented across five human mTEC samples, with mature mTECs showing increased rates of global transcript mis-initiation. Both mTEC populations have increased splicing entropy, which appears to be driven by expression of peripheral splicing factors. Furthermore, EREs enriched in long terminal repeat retrotransposons are up-regulated during mTEC maturation and enriched in genomic proximity to differentially expressed genes. We provide an interactive interface to explore the transcriptome diversity we uncovered at http://transcriptomediversity.cshl.edu/. Our findings represent an important first step towards the generation of a comprehensive map of transcriptome diversity in the healthy human thymus. Ultimately, a complete map of thymic expression diversity will allow for the identification of epitopes that contribute to the pathogenesis of auto-immunity and that drive immune recognition of tumor antigens.


## This repository
This repository contains all analyses conducted in the manuscript.
It is divided into subdirectories(brief description below), where
each subdirectory contains analysis organised in either a snakemake workflow
(with the workflow depicted as a DAG in the respective README.md) or
consecutively labelled scripts. 

All manuscript figures can be reproduced by running the python notebooks in the 
Figures folder. To make this possible, source data for all the figures are provided,
some of which exceed the standard github file size. To ensure proper file download,
please use [git lsf](https://git-lfs.github.com/) for cloning this repository. 


### 5P-Cap-seq-alignment
* Alignment and quality control of 5P-Cap-seq data of mTEC samples from 5 human tissue donors

### 5P-Cap-seq-processing
* Transcription start site and transcription start region calling for immature and mature human mTECs

### Gene-lists
* gene lists derived from public datasets, and re-mapping to GRCh38
* estimating tissue specifc antigens based on tissue exclusivity index tau and
    22 tissues from GTEx
    
### ERE-analysis
ERE expression in medullary thymic epithelial cells.

## Dependencies
- Python: python (v3.9.6), IPython (v7.26.0), scipy (v1.7.1), seaborn (v0.11.1), matplotlib (v3.4.2), matplotlib_venn (v0.11.6), numpy (v1.21.1), pandas (v1.3.1),  sklearn (v0.24.2), statsmodels (v0.12.2)
- R: R (v4.0.3), CAGEr (v1.32), limma (v3.46), Sleuth (v0.30.0), biomaRt (v2.46.3), rMATS (v4.1.1)
- sequence analysis: umi_tools (v1.1), fastq_screen (v0.14.0), STAR (v2.7.2b), samtools (v1.11), picard (v2.18.20), multiqc (v1.9), paraclu (v9), bedops (v2.4.38), bedtools (v2.29.2), HOMER (v4.11.1), sambamba (v0.8), fastp (v0.11.8), Kallisto (v0.4.6), ngs.plot.r (v2.61), clustal omega (v1.2.4)
