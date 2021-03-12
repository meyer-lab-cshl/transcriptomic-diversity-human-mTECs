#################################################################
# Differential expression with DESeq2
#################################################################

library(dplyr)
library(tidyr)
library(DESeq2)

setwd('/Users/mpeacey/TE_thymus/analysis/hi_vs_lo_local')

## Set parameters

p_value_cutoff = 0.05
log_fold_change_cutoff = 0.58

## Import count matrix

data = read.table("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo_local/count_table/TE_local_hi_vs_lo.cntTable",
                  header=T,row.names=1)
colnames(data) = c('214_HI', '214_LO', '221_HI', '221_LO', '226_HI', '226_LO')

## Subset into TE count matrix 'TE_data'

TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

## Filter count matrix to exclude non-expressed genes

min_read = 1
TE_data = TE_data[apply(TE_data,1,function(x){max(x)}) > min_read,]

## Define sampleInfo

ID = colnames(TE_data)
sampleInfo = data.frame(ID,row.names=colnames(TE_data))
sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'group'), sep = '_'))
sampleInfo$patient = factor(sampleInfo$patient)

## Construct DESeq dataset object

dds <- DESeqDataSetFromMatrix(countData = TE_data, colData = sampleInfo, design = ~patient + group)
dds$group = relevel(dds$group,ref="LO")

## Run differential expression analysis

dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)

results_df = as.data.frame(res)
