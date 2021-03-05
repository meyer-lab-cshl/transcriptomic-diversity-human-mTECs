library(dplyr)
library(readr)
library(tidyr)
library(DESeq2)

## Import count matrix

data = read.table("hi_vs_lo.cntTable",header=T,row.names=1)

## Subset into 'gene' and 'TE' count matrices

gene_data = data[grep("^ENSG",rownames(data)),]
TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

## Convert 'gene_ID' to 'gene_name' from GTF file

GENCODE_annotation = read.table(file = 'gencode.v38_gene_annotation_table.txt', header = 1)

## Filter count matrix to exclude non-expressed genes

min_read = 1
gene_data <- gene_data[apply(gene_data,1,function(x){max(x)}) > min_read,]

## Define sampleInfo

ID = colnames(gene_data)
sampleInfo = data.frame(ID,row.names=colnames(gene_data))
sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'group'), sep = '_'))

## Modify gene_data

gene_data <- cbind(Geneid = rownames(gene_data), gene_data)
gene_data = merge(gene_data, GENCODE_annotation, by = 'Geneid')
#rownames(gene_data) = gene_data$GeneSymbol 

# Doesn't work because some gene IDs are not unique!

## Construct DESeq dataset object

dds <- DESeqDataSetFromMatrix(countData = gene_data, colData = sampleInfo, design = ~group+patient)
#dds$groups = relevel(dds$groups,ref="mTEC_lo")

## Run differential expression analysis

dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)

## Transform raw count data

vs_dds <- vst(dds, blind=FALSE)

## Export as csv

write.table(res, file="hi_vs_lo_gene_TE_analysis.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file="hi_vs_lo_sigdiff_gene_TE.txt",sep="\t", quote=F)
