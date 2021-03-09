library(dplyr)
library(readr)
library(tidyr)
library(DESeq2)
library(reshape2)
library(ggplot2)
library(svglite)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(pheatmap)

## Set parameters

p_value_cutoff = 0.05
log_fold_change_cutoff = 1

#################################################################
# Differential expression with DESeq2
#################################################################

## Import count matrix

data = read.table("hi_vs_lo.cntTable",header=T,row.names=1)
colnames(data) = c('214_HI', '221_HI', '226_HI', '214_LO', '221_LO', '226_LO')

## Subset into TE count matrix 'TE_data'

TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

## Filter count matrix to exclude non-expressed genes

min_read = 1
TE_data = TE_data[apply(TE_data,1,function(x){max(x)}) > min_read,]

## Define sampleInfo

ID = colnames(TE_data)
sampleInfo = data.frame(ID,row.names=colnames(TE_data))
sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'group'), sep = '_'))

## Construct DESeq dataset object

dds <- DESeqDataSetFromMatrix(countData = TE_data, colData = sampleInfo, design = ~patient + group)
dds$group = relevel(dds$group,ref="LO")

## Run differential expression analysis

dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)

## Convert results to dataframe and add signficance label

results_df = as.data.frame(res)

results_df = cbind(ID = rownames(results_df), results_df)
results_df = separate(data = results_df, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
results_df = cbind(ID = rownames(results_df), results_df)

results_df = mutate(results_df, significant = case_when((padj < p_value_cutoff) & (abs(log2FoldChange) > log_fold_change_cutoff) ~ TRUE, 
                                        (padj >= p_value_cutoff) | (abs(log2FoldChange) <= log_fold_change_cutoff) ~ FALSE))

results_df = mutate(results_df, abs_log2FoldChange = abs(log2FoldChange))

sig_results_df = results_df[results_df$significant == TRUE,]

## Transform raw count data and convert to dataframe

vs_dds <- vst(dds, blind=FALSE)
transformed_counts = as.data.frame(assay(vs_dds))

sigGenes = rownames(results_df[results_df$significant == TRUE,])
sig_transformed_counts = assay(vs_dds)[rownames(normalized_counts) %in% sigGenes,]

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

PCA = plotPCA(vs_dds, intgroup = 'group') + 
  ggtitle('mTEC-hi vs mTEC-lo', 'TE expression PCA') +
  geom_text(aes(label = colnames(vs_dds)), nudge_x = 0.5, nudge_y = 0.2)

PCA + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                        plot.subtitle = element_text(size = 14),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14),
                        axis.title = element_text(size = 14),
                        axis.line = element_line(size = 0.8),
                        panel.border = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_PCA.png", 
       width = 20, height = 15, units = "cm")

#################################################################
# Heatmap version 1 (ggplot)
#################################################################

sig_normalized_counts_long <- melt(sig_normalized_counts, id.vars=c("ID"))

ggplot(sig_normalized_counts_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_raster() + theme(axis.text.x=element_text(angle=65, hjust=1)) 

#################################################################
# Heatmap version 2 (pheatmap)
#################################################################

## Rows ordered by fold change w/out clustering

select <- order(sig_results_df$abs_log2FoldChange, decreasing=TRUE)
pheatmap(sig_transformed_counts[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)

## Rows clustered

my_heatmap = pheatmap(sig_transformed_counts, 
                      cluster_rows=TRUE,
                      show_rownames=TRUE, 
                      cluster_cols=TRUE,
                      cutree_cols = 2)

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_heatmap.png")

#################################################################
# Volcano
#################################################################

volcano_plot = ggplot(data = results_df, aes(x = log2FoldChange, y = -log10(padj), colour = significant)) +
  geom_point(alpha = 0.6, aes(colour = significant)) +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  xlim(-3, 3) +
  ggtitle('mTEC-hi vs mTEC-lo', 'TE expression') +
  scale_colour_manual(values = c('#9B9A99', "red")) +
  guides(colour = FALSE) +
  geom_label_repel(
    data = subset(results_df, significant == TRUE),
    aes(label = subset(results_df, significant == TRUE)$gene),
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

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_volcano_plot.png", 
       width = 20, height = 15, units = "cm")

#################################################################
# Pie chart
#################################################################

## Frequency of all TEs in mTEC-HI cells

normalized_counts_HI = as.data.frame(counts(dds, normalized = TRUE)) %>%
select(c('214_HI', '221_HI', '226_HI'))
normalized_counts_HI$mean = rowMeans(normalized_counts_HI)
normalized_counts_HI = select(normalized_counts_HI, mean)

normalized_counts_HI = cbind(ID = rownames(normalized_counts_HI), normalized_counts_HI)
normalized_counts_HI = separate(data = normalized_counts_HI, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
normalized_counts_HI = cbind(ID = rownames(normalized_counts_HI), normalized_counts_HI)

normalized_counts_HI = mutate(normalized_counts_HI, class = sub("\\?", "", class))

normalized_counts_HI_class = group_by(normalized_counts_HI, class) %>% summarize(sum = sum(mean))
normalized_counts_HI_class = as.data.frame(normalized_counts_HI_class)

ggplot(normalized_counts_HI_class, aes(x="", y=sum, fill=class)) + 
  geom_col(poisiton = 'fill')


######
sig_normalized_counts = as.data.frame(sig_normalized_counts)
sig_normalized_counts_HI = select(sig_normalized_counts, c('214_HI', '221_HI', '226_HI'))
sig_normalized_counts_HI$mean = rowMeans(sig_normalized_counts_HI)
sig_normalized_counts_HI = select(sig_normalized_counts_HI, mean)

sig_normalized_counts_HI = cbind(ID = rownames(sig_normalized_counts_HI), sig_normalized_counts_HI)
sig_normalized_counts_HI = separate(data = sig_normalized_counts_HI, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
sig_normalized_counts_HI = cbind(ID = rownames(sig_normalized_counts_HI), sig_normalized_counts_HI)

sig_normalized_counts_HI = merge(sig_normalized_counts_HI, 
                                 sig_results_df)

sig_normalized_counts_HI = select(sig_normalized_counts_HI, c('gene', 'family', 'class', 'mean'))

sig_normalized_counts_HI_class = group_by(sig_normalized_counts_HI, class) %>% summarize(sum = sum(mean))
sig_normalized_counts_HI_class = as.data.frame(sig_normalized_counts_HI_class)

ggplot(sig_normalized_counts_HI_class, aes(x="", y=sum, fill=class)) + 
  geom_bar(stat="identity", width=1)




