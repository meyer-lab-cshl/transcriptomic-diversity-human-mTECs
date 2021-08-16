library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

#################################################################
# PCA w/ GTEx data (A)
#################################################################

pcaData = plotPCA(vs_dds_transcripts_TE, intgroup='tissue', returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

PCA = ggplot(pcaData, aes(PC1, PC2, fill = tissue)) + 
  geom_point(size=3, shape = 21, stroke = 0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  labs(fill= "Tissue") +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff')) 

PCA + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                         plot.subtitle = element_text(size = 14),
                         axis.text.x = element_text(size = 11),
                         axis.text.y = element_text(size = 11),
                         axis.title.x = element_text(size = 13, margin = margin(t = 8)),
                         axis.title.y = element_text(size = 13),
                         axis.line = element_line(size = 0.8),
                         panel.border = element_blank(),
                         legend.text = element_text(size = 10),
                         legend.title = element_text(size = 13))

#################################################################
# Volcano plot (B)
#################################################################

input = results_df_transcripts_TE

input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))

input = mutate(input, grouped_class = case_when(class == 'DNA' ~ 'DNA', 
                                                 class == 'LINE' ~ 'LINE',
                                                 class == 'LTR' ~ 'LTR',
                                                 class == 'SINE' ~ 'SINE',
                                                 class == 'Retroposon' ~ 'Other',
                                                 class == 'Satellite' ~ 'Other'))

## LTRs highlighted

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = grouped_class), size = 2.5, alpha = 0.8, shape = 21, stroke = 0.4) +
  geom_point(data = subset(input, significant == FALSE), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
  xlab(expression('log'[2]*' fold-change')) +
  ylab(expression('-log'[10]*' p-value')) +
  xlim(-3, 3) +
  scale_fill_brewer(palette = 'Set1') +
  labs(fill= "Class")

#ggtitle('mTEC-hi vs mTEC-lo', 'TE transcripts') 

## Plot

volcano_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank(),
                                  legend.text = element_text(size = 11),
                                  legend.title = element_text(size = 14),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/transcripts_volcano.png", 
       width = 11, height = 5.25, units = "in")
