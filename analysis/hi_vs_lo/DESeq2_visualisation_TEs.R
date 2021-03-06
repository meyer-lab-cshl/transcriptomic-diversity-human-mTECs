library(ggplot2)
library(svglite)
library(ggpubr)
library(ggrepel)

##This is a set of plot templates for visualizing DESeq2 output. It assumes
##the following assignments:
##
## DESeq dataset object: dds
## Differential expression results: res
## Transformed counts: vs_dds

## PCA

plotPCA(vs_dds, intgroup = 'group') + theme_bw()

## Count plots for individual genes

par(mfrow=c(2,2))
plotCounts(dds, gene = 'ENSG00000160224.17', intgroup = 'group') ## AIRE
plotCounts(dds, gene = 'ENSG00000121594.12', intgroup = 'group') ## CD80
plotCounts(dds, gene = 'ENSG00000153266.13', intgroup = 'group') ## FEZF2
plotCounts(dds, gene = '', intgroup = 'group') ## FEZF2

## Volcano

df2 = filter(df, Class == 'protein_coding')

library(EnhancedVolcano)

EnhancedVolcano(df,
                lab = rownames(df),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'mTEC-lo vs mTEC-hi gene expression',
                pCutoff = 0.05)

volcano_plot = ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), colour = significant)) +
                geom_point(alpha = 0.6) +
                geom_vline(xintercept = c(1, -1), linetype = 'dashed') +
                geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
                xlab(expression('Log'[2])) +
                ylab(expression('-Log'[10])) +
                xlim(-3, 3) +
                ggtitle('mTEC-lo vs mTEC-hi', 'TE expression') +
                scale_colour_manual(values = c('#9B9A99', "red")) +
                guides(colour = FALSE) +
                geom_text_repel(
                data = subset(df, significant == TRUE),
                aes(label = rownames(subset(df, significant == TRUE))),
                size = 4)

volcano_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/lo_vs_hi_TEs_volcano_plot.png", 
       width = 20, height = 15, units = "cm")
