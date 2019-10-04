---
title: "Gviz script for genomic ranges"
author: "Léonie Strömich"
date: "August 31, 2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Script to calculate distribution of tags over gene body. Each TSS is counted to a certain gene body feature class exactely once. 

1. Define overlap of TSS with gene body genomic ranges. Each TSS is counted exactely once depending on the following priorities:
    a. TSS_10 
    b. TSS_100
    c. first_exon 
    d. TSS_10_other_exons
    e. TSS_100_other_exons
    f. TTS_100
    g. upstream
    h. other_exon
    i. intron
    j. antisense
    k. downstream

2. Add count of TSS to the gene class, the gene of the found feature range, belongs to. 

These steps are performed for TSS in mTECs and tissues seperately. The resulting matrices are merged later and used for plotting TSS distribution in different gene classes. 

Input: Gene body genomic ranges, Gene classes, Filtered positions in mTECs and tissues

Output: Tag count matrices with gene classes as rows and Positions as columns

Running time: ~6h for filtered datasets, ~7h for unfiltered dataset

## Prerequisites
#### Load libraries

```{r libraries, message=FALSE}
library(Gviz)
```

#### Functions
```{r functions, eval=TRUE}
preprocess_df <- function(x) {
  x[x$strand == "-",]$BarcodeCount <- (x[x$strand == "-",]$BarcodeCount)*-1
  x$strand <- "+"
  return(x)
}

```

#### Read in data 

```{r read in, eval = FALSE}
### read in  --> Transfer/tss_data/Filtered
mTEC_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered.csv", header=TRUE, sep=",")
tissues_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.filtered.csv", header=TRUE, sep=",")

```

## Analysis

Pre-processing
```{r pre_processing, eval=FALSE}
mTEC_pos_fil <- preprocess_df(mTEC_pos_fil)
mTEC_ranges <- makeGRangesFromDataFrame(mTEC_pos_fil, keep.extra.columns = TRUE)
tissues_pos_fil <- preprocess_df(tissues_pos_fil)
tissues_ranges <- makeGRangesFromDataFrame(tissues_pos_fil, keep.extra.columns = TRUE)


## list of genes, corresponding chromosomes and genomic coordinates 
genes <- c("CTRB2", "CYP3A4","USHBP1","STAR","LDHC","CELA2B","RARRES1","CHIT1","RHCG","CD46","TG")
chromosomes <- list("CTRB2" = "chr16", "CYP3A4"="chr7","USHBP1"="chr19","STAR"="chr8","LDHC"="chr11",
                    "CELA2B"="chr1","RARRES1"="chr3","CHIT1"="chr1","RHCG"="chr15","CD46"="chr1","TG"="chr8")
                    
starts <- list("CTRB2" = 75237550, "CYP3A4"=99351650,"USHBP1"=17359350,"STAR"=37998880,"LDHC"=18430380,
               "CELA2B"=15800650,"RARRES1"=158412450,"CHIT1"=203184000,"RHCG"=90012230,"CD46"=207921900,"TG"=133862210)
ends <- list("CTRB2" = 75241200, "CYP3A4"=99382900,"USHBP1"=17376000,"STAR"=38008780,"LDHC"=18473310,
             "CELA2B"=15819850,"RARRES1"=158451900,"CHIT1"=203199350,"RHCG"=90040450,"CD46"=207969400,"TG"=134150850)

```

Loop over genes in list and plot genomic ranges. 

```{r loop, eval=TRUE}
for(i in genes[3]){ # change here 
  #read in genomic location from list
  chromosome <- chromosomes[[i]]
  afrom <- starts[[i]]
  ato <- ends[[i]]
  
  
  m_max <- max(mTEC_pos_fil[mTEC_pos_fil$chrom == chromosome & mTEC_pos_fil$start > afrom &  mTEC_pos_fil$start < ato ,]$BarcodeCount)
  if(m_max < 0) {m_max = 0}
  m_min <- min(mTEC_pos_fil[mTEC_pos_fil$chrom == chromosome & mTEC_pos_fil$start > afrom &  mTEC_pos_fil$start < ato ,]$BarcodeCount)
  if(m_min > 0) {m_min = 0}
  mtrack <- DataTrack(mTEC_ranges, type="h", baseline=0, name="mTECs", 
                      col.baseline="grey", ylim=c(m_min-(0.1*m_min),m_max+(0.1*m_max)))
  
  t_max <- max(tissues_pos_fil[tissues_pos_fil$chrom == chromosome & tissues_pos_fil$start > afrom &  tissues_pos_fil$start < ato ,]$BarcodeCount)
  if(t_max < 0) {t_max = 0}
  t_min <- min(tissues_pos_fil[tissues_pos_fil$chrom == chromosome & tissues_pos_fil$start > afrom &  tissues_pos_fil$start < ato ,]$BarcodeCount)
  if(t_min > 0) {t_min = 0}
  ttrack <- DataTrack(tissues_ranges, type="h", baseline=0, name="tissues wo thymus", 
                      col.baseline="grey", ylim=c(t_min-(0.1*t_min),t_max+(0.1*t_max)))
  
  ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chromosome)
  axisTrack <- GenomeAxisTrack()
  
  
  ## Create genetrack from UCSC
  # fill = "#8282d2" for colours 
  for(x in c("dense","squish")){
    genetrack=UcscTrack(track="RefSeq Genes", table="refGene", trackType="GeneRegionTrack",chromosome=chromosome, genome="hg19",
                      stacking = x,
                      rstart="exonStarts", rends="exonEnds", gene="name", symbol="name2", transcript="name", strand="strand", 
                      name="RefSeq Genes", feature="name2", geneSymbol=T, from = afrom, to = ato)
  
  
  
  plot_file_name = paste("/home/stroemic/hiwi_16/analysis/shifted_TSS/plots/gviz/",i,"_gviz_range_",x,"_.pdf",sep="")
  #pdf(plot_file_name, width=12)
  plotTracks(list(ideoTrack, axisTrack, genetrack, ttrack, mtrack), from = afrom, to = ato,
             showBandId = TRUE, showID=TRUE, showFeatureID=TRUE,transcriptAnnotation = "symbol"
             groupAnnotation = "id")
  #dev.off()
  }
  
  
}

```


