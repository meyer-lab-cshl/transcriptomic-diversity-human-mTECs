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
log_fold_change_cutoff = 0.58

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

results_df = mutate(results_df, significant = case_when(padj < p_value_cutoff ~ TRUE, padj >= p_value_cutoff ~ FALSE))

results_df = mutate(results_df, FC_significant = case_when(abs(log2FoldChange) > log_fold_change_cutoff ~ TRUE, 
                                                           abs(log2FoldChange) <= log_fold_change_cutoff ~ FALSE))

results_df = mutate(results_df, abs_log2FoldChange = abs(log2FoldChange))

sig_results_df = results_df[results_df$significant == TRUE,]

## Transform raw count data 

vs_dds <- vst(dds, blind=FALSE)
transformed_counts = as.data.frame(assay(vs_dds))

sigGenes = rownames(results_df[results_df$significant == TRUE,])
sig_transformed_counts = assay(vs_dds)[rownames(transformed_counts) %in% sigGenes,]

upGenes = rownames(results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange > 0),])

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
# Stacked bar chart: class/family frequency
#################################################################

build_count_table = function(group, mode){
  
  normalized_counts = as.data.frame(counts(dds, normalized = TRUE))
  
  for (i in group){
    
    if (i == 'all'){
      
      input = normalized_counts
      
    }
    
    else if (i == 'differentially_regulated'){
      
      input = normalized_counts[rownames(normalized_counts) %in% sigGenes,]
      
    }
    
    else if (i == 'upregulated'){
      
      input = normalized_counts[rownames(normalized_counts) %in% upGenes,]
      
    }
    
    input_HI = select(input, c('214_HI', '221_HI', '226_HI'))
    input_HI$mean = rowMeans(input_HI)
    input_HI = select(input_HI, mean)
    
    input_HI = cbind(ID = rownames(input_HI), input_HI)
    input_HI = separate(data = input_HI, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
    input_HI = cbind(ID = rownames(input_HI), input_HI)
    
    input_HI = mutate(input_HI, class = sub("\\?", "", class))
    
    if (mode == 'class'){
      
      input_HI_class = group_by(input_HI, class) %>% summarize(sum = sum(mean))
      input_HI_class = as.data.frame(input_HI_class)
      input_HI_class$group = i
      
      output = input_HI_class
      
    }
    
  }

  return(output)
      
}

normalized_counts_HI_LTR_family = normalized_counts_HI %>% 
                                  filter(class == 'LTR') %>%
                                  group_by(family) %>% 
                                  summarize(sum = sum(mean))
normalized_counts_HI_LTR_family = as.data.frame(normalized_counts_HI_LTR_family)
normalized_counts_HI_LTR_family$group = 'All'

## Frequency of differentially expressed TEs in mTEC-HI cells

normalized_counts = as.data.frame(counts(dds, normalized = TRUE))
sig_normalized_counts = normalized_counts[rownames(normalized_counts) %in% sigGenes,]

sig_normalized_counts_HI = select(sig_normalized_counts, c('214_HI', '221_HI', '226_HI'))
sig_normalized_counts_HI$mean = rowMeans(sig_normalized_counts_HI)
sig_normalized_counts_HI = select(sig_normalized_counts_HI, mean)

sig_normalized_counts_HI = cbind(ID = rownames(sig_normalized_counts_HI), sig_normalized_counts_HI)
sig_normalized_counts_HI = separate(data = sig_normalized_counts_HI, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
sig_normalized_counts_HI = cbind(ID = rownames(sig_normalized_counts_HI), sig_normalized_counts_HI)

sig_normalized_counts_HI = mutate(sig_normalized_counts_HI, class = sub("\\?", "", class))

sig_normalized_counts_HI_class = group_by(sig_normalized_counts_HI, class) %>% summarize(sum = sum(mean))
sig_normalized_counts_HI_class = as.data.frame(sig_normalized_counts_HI_class)
sig_normalized_counts_HI_class$group = 'Differentially regulated'

sig_normalized_counts_HI_LTR_family = sig_normalized_counts_HI %>% 
                                      filter(class == 'LTR') %>%
                                      group_by(family) %>% 
                                      summarize(sum = sum(mean))

sig_normalized_counts_HI_LTR_family = as.data.frame(sig_normalized_counts_HI_LTR_family)
sig_normalized_counts_HI_LTR_family$group = 'Differentially regulated'

## Frequency of upregulated TEs in mTEC-HI cells

normalized_counts = as.data.frame(counts(dds, normalized = TRUE))
up_normalized_counts = normalized_counts[rownames(normalized_counts) %in% upGenes,]

up_normalized_counts_HI = select(up_normalized_counts, c('214_HI', '221_HI', '226_HI'))
up_normalized_counts_HI$mean = rowMeans(up_normalized_counts_HI)
up_normalized_counts_HI = select(up_normalized_counts_HI, mean)

up_normalized_counts_HI = cbind(ID = rownames(up_normalized_counts_HI), up_normalized_counts_HI)
up_normalized_counts_HI = separate(data = up_normalized_counts_HI, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
up_normalized_counts_HI = cbind(ID = rownames(up_normalized_counts_HI), up_normalized_counts_HI)

up_normalized_counts_HI = mutate(up_normalized_counts_HI, class = sub("\\?", "", class))

up_normalized_counts_HI_class = group_by(up_normalized_counts_HI, class) %>% summarize(sum = sum(mean))
up_normalized_counts_HI_class = as.data.frame(up_normalized_counts_HI_class)
up_normalized_counts_HI_class$group = 'Upregulated'

up_normalized_counts_HI_LTR_family = up_normalized_counts_HI %>% 
  filter(class == 'LTR') %>%
  group_by(family) %>% 
  summarize(sum = sum(mean))

up_normalized_counts_HI_LTR_family = as.data.frame(up_normalized_counts_HI_LTR_family)
up_normalized_counts_HI_LTR_family$group = 'Upregulated'

## Append the dataframes

class_frequencies = bind_rows(normalized_counts_HI_class, 
                              sig_normalized_counts_HI_class, 
                              up_normalized_counts_HI_class)
LTR_family_frequencies = bind_rows(normalized_counts_HI_LTR_family, 
                                   sig_normalized_counts_HI_LTR_family,
                                   up_normalized_counts_HI_LTR_family)

## Plot stacked bars

bar_chart = ggplot(class_frequencies, aes(x = group, y = sum, fill = class)) + 
            geom_col(colour = 'black', position = 'fill') +
            scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
            scale_fill_brewer(palette = "Set1") +
            xlab('') +
            ylab('Fraction of normalized reads')

bar_chart + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                               plot.subtitle = element_text(size = 14),
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_line(color = 'grey'),
                               panel.grid.minor.y = element_blank(),
                               axis.text.x = element_text(size = 14),
                               axis.text.y = element_text(size = 14),
                               axis.title = element_text(size = 14),
                               axis.line = element_line(size = 0.8),
                               panel.border = element_blank(),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 14))

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_classbreakdown.png", 
       width = 15, height = 15, units = "cm")