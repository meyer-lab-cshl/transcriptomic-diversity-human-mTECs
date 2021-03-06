library(dplyr)
library(readr)
library(tidyr)
library(DESeq2)

## Import count matrix

data = read.table("hi_vs_lo.cntTable",header=T,row.names=1)
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
dds$group = relevel(dds$group,ref="lo")

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

df = mutate(df, significant = case_when((padj < 0.05) & (abs(log2FoldChange) > 1) ~ TRUE, 
                                        (padj >= 0.05) | (abs(log2FoldChange) <= 1) ~ FALSE))

## Export

#write.table(res, file="hi_vs_lo_gene_TE_analysis.txt", sep="\t",quote=F)
#resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
#write.table(resSig, file="hi_vs_lo_sigdiff_gene_TE.txt",sep="\t", quote=F)
