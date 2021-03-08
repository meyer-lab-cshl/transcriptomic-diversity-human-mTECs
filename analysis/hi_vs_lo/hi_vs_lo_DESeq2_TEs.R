library(dplyr)
library(readr)
library(tidyr)
library(DESeq2)
library(reshape2)
library(ggplot2)
library(svglite)
library(ggpubr)
library(ggrepel)
library(viridis)
library(pheatmap)

#################################################################
# Differential expression with DESeq2
#################################################################

## Import count matrix

data = read.table("hi_vs_lo.cntTable",header=T,row.names=1)

## Subset into TE count matrix 'TE_data'

TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

## Filter count matrix to exclude non-expressed genes

min_read = 1
TE_data <- TE_data[apply(TE_data,1,function(x){max(x)}) > min_read,]

## Define sampleInfo

ID = colnames(TE_data)
sampleInfo = data.frame(ID,row.names=colnames(TE_data))
sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'group'), sep = '_'))

## Construct DESeq dataset object

dds <- DESeqDataSetFromMatrix(countData = TE_data, colData = sampleInfo, design = ~patient + group)
dds$group = relevel(dds$group,ref="lo")

## Run differential expression analysis

dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)

## Convert results to dataframe and add signficance label

results_df = as.data.frame(res)

results_df = cbind(ID = rownames(results_df), results_df)
results_df = separate(data = results_df, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
results_ID = cbind(ID = rownames(results_df), results_df)

results_df = mutate(results_df, significant = case_when((padj < 0.05) & (abs(log2FoldChange) > 1) ~ TRUE, 
                                        (padj >= 0.05) | (abs(log2FoldChange) <= 1) ~ FALSE))

results_df = mutate(results_df, abs_log2FoldChange = abs(log2FoldChange))

sig_results_df = results_df[results_df$significant == TRUE,]

## Transform raw count data and convert to dataframe

vs_dds <- vst(dds, blind=FALSE)
normalized_counts = as.data.frame(assay(vs_dds))

sigGenes = rownames(results_df[results_df$significant == TRUE,])
sig_normalized_counts = assay(vs_dds)[rownames(normalized_counts) %in% sigGenes,]

## Export

#write.table(res, file="hi_vs_lo_gene_TE_analysis.txt", sep="\t",quote=F)
#resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
#write.table(resSig, file="hi_vs_lo_sigdiff_gene_TE.txt",sep="\t", quote=F)

#################################################################
#################################################################
# Data visualization
#################################################################
#################################################################

#################################################################
# PCA
#################################################################

plotPCA(vs_dds, intgroup = 'group') + theme_bw()

#################################################################
# Heatmap version 1
#################################################################

normalized_counts = as.data.frame(assay(vs_dds))
normalized_counts = cbind(ID = rownames(normalized_counts), normalized_counts)
#normalized_counts = separate(data = normalized_counts, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')

sigGenes = rownames(df[df$significant == TRUE,])
sig_normalized_counts <- normalized_counts[rownames(normalized_counts) %in% sigGenes,]

sig_normalized_counts_long <- melt(sig_normalized_counts, id.vars=c("ID"))

ggplot(sig_normalized_counts_long, aes(x = variable, y = ID, fill = value)) +
  geom_raster() + theme(axis.text.x=element_text(angle=65, hjust=1)) +
  scale_fill_viridis(trans = 'sqrt')

#################################################################
# Heatmap version 2
#################################################################

## Rows ordered by fold change w/out clustering

select <- order(sig_results_df$abs_log2FoldChange, decreasing=TRUE)
pheatmap(sig_normalized_counts[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)

## Rows clustered

pheatmap(sig_normalized_counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE)


#################################################################
# Volcano
#################################################################

volcano_plot = ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), colour = significant)) +
  geom_point(alpha = 0.6, aes(colour = significant)) +
  xlab(expression('Log'[2]*' FC (HI/LO)')) +
  ylab(expression('-Log'[10]*' P value')) +
  xlim(-3, 3) +
  ggtitle('mTEC-lo vs mTEC-hi', 'TE expression') +
  scale_colour_manual(values = c('#9B9A99', "red")) +
  guides(colour = FALSE) +
  geom_label_repel(
    data = subset(df, significant == TRUE),
    aes(label = subset(df, significant == TRUE)$class),
    size = 3.5, 
    box.padding = unit(0.9, 'lines'), 
    point.padding = unit(0.2, 'lines'),
    max.overlaps = 20)

#geom_vline(xintercept = c(1, -1), linetype = 'dashed') +
#geom_hline(yintercept = -log10(0.05), linetype = 'dashed')

volcano_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/lo_vs_hi_TEs_volcano_plot.png", 
       width = 20, height = 15, units = "cm")

