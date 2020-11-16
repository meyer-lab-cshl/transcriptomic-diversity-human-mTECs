# Genome-wide TSS in medullary TEC and peripheral tissues

Antigen epitopes missing in the thymus can result in autoimmunity by allowing auto-reactive T cell escape into the periphery.
Here, we will generate comprehensive epitope maps of the human thymus. Using 5Pseq and RNAseq, we will first map transcription
start sites and splicing events in the thymus and compare them to corresponding maps for peripheral tissues (FANTOM data).
This comparison will allow us to detect mis-initiation and splicing events leading to potentially altered epitope maps in the thymus.

This repository contains the data processing steps, from RNAseq reads to transcription start site maps. It is divided by data
into subdirectories. Each subdirectory contains analysis organised in a snakemake workflow, with the workflow depicted as a DAG in the 
respective README.md

## human
* Analysis of 5Pseq data from 5 human tissue donors

## mouse
* Analysis of 5Pseq data from mouse embryonic stem cells

## fantom
* Processing of human and mouse fantom5 samples
