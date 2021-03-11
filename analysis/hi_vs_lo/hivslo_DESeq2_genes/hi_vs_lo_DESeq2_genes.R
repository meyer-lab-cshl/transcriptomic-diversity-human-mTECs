#################################################################
# Differential expression with DESeq2
#################################################################

library(dplyr)
library(readr)
library(tidyr)
library(DESeq2)
library(reshape2)
library(svglite)
library(gridExtra)
library(pheatmap)

## Set parameters

p_value_cutoff = 0.05
log_fold_change_cutoff = 0.58

## Import count matrix

setwd("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo")
data = read.table("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/hi_vs_lo.cntTable",header=T,row.names=1)
colnames(data) = c('214_HI', '221_HI', '226_HI', '214_LO', '221_LO', '226_LO')

## Subset into gene count matrix 'gene_data'

gene_data = data[grep("^ENSG",rownames(data)),]

## Filter count matrix to exclude non-expressed genes

min_read = 1
gene_data <- gene_data[apply(gene_data,1,function(x){max(x)}) > min_read,]

## Define sampleInfo

ID = colnames(gene_data)
sampleInfo = data.frame(ID,row.names=colnames(gene_data))
sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'group'), sep = '_'))

## Construct DESeq dataset object

dds <- DESeqDataSetFromMatrix(countData = gene_data, colData = sampleInfo, design = ~patient + group)
dds$group = relevel(dds$group,ref="LO")

## Run differential expression analysis

dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)

## Transform raw count data

vs_dds <- vst(dds, blind=FALSE)

## Add GTF annotation information to results table

df = as.data.frame(res)

GENCODE_annotation = read.table(file = 'gencode.v38_gene_annotation_table.txt', header = 1)
GENCODE_annotation_subset = select(GENCODE_annotation, -c(Start, End, Strand, Length))

df = cbind(Geneid = rownames(df), df)
df = merge(df, GENCODE_annotation_subset, by = 'Geneid')

## Add p value annotation

df = mutate(df, significant = case_when(padj < p_value_cutoff ~ TRUE, padj >= p_value_cutoff ~ FALSE))

df = mutate(df, FC_significant = case_when(abs(log2FoldChange) > log_fold_change_cutoff ~ TRUE, 
                                                           abs(log2FoldChange) <= log_fold_change_cutoff ~ FALSE))

df = mutate(df, overall_significant = case_when((significant == TRUE) & (FC_significant == TRUE) ~ TRUE, 
                                                                (significant == FALSE) | (FC_significant == FALSE) ~ FALSE))