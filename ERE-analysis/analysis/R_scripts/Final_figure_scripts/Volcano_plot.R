library(tidyverse)

working_directory = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis'

################################################################################
# TE transcripts
################################################################################

input = readRDS(file = glue('{working_directory}/R_variables/results_df_transcripts_ERE')) %>%
  


input$class = factor(input$class, levels = c('LTR', 'SINE', 'LINE', 'Retroposon'))

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6), stroke = 0) +
  geom_point(data = subset(input, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = class), size = 1.8, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, significant == FALSE), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  xlim(-2, 3) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#F781BF')) +
  labs(fill= "")

volcano_plot = ggplot() +
  geom_point(data = subset(input, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = class), size = 1.8, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, significant == FALSE), aes(x = log2FoldChange, y = -log10(padj)), size = 1.8, shape = 16, stroke = 0, color = alpha('#9B9A99', 0.6)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  xlim(-2, 3) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#F781BF')) +
  labs(fill= "")

## ereMAPs annotated

input = mutate(input, ereMAP = case_when(gene %in% ereMAPs$ERE.family ~ T,
                                         !(gene %in% ereMAPs$ERE.family) ~ F))

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, ereMAP == T & significant == T), aes(x = log2FoldChange, y = -log10(padj), fill = significant), size = 2, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, ereMAP == F), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#fb9a99')) +
  labs(fill= "") +
  ggrepel::geom_label_repel(data = filter(input, ereMAP == T & significant == T),
                            aes(x = log2FoldChange,
                                y = -log10(padj), 
                                label = gene),
                            force = 2)

## LIONS annotated

input = mutate(input, LIONS = case_when(gene %in% subset(LIONS, total_occurences >= 2)$`sub-family` ~ T,
                                        !(gene %in% subset(LIONS, total_occurences >= 2)$`sub-family`) ~ F)) %>%
  

#input = mutate(input, LIONS = case_when(gene %in% subset(LIONS, classification == 'HI')$`sub-family` ~ 'HI', 
#                                        gene %in% subset(LIONS, classification == 'LO')$`sub-family` ~ 'LO',
#                                        !(gene %in% LIONS$`sub-family`) ~ 'NONE'))

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, LIONS == T & significant == T), aes(x = log2FoldChange, y = -log10(padj), fill = LIONS), size = 2, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, LIONS == F), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#fb9a99')) +
  labs(fill= "") +
  ggrepel::geom_label_repel(data = filter(input, LIONS == T), 
                           aes(x = log2FoldChange, 
                               y = -log10(padj), 
                               label = gene))

#################################################################
# TE local
#################################################################

input = results_df_local_ERE

input$class = factor(input$class, levels = c('LTR', 'SINE', 'LINE', 'Retroposon'))

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = class), size = 2, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, significant == FALSE), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#fb9a99')) +
  labs(fill= "")

## ereMAPs annotated

#input = mutate(input, ereMAP = case_when(locus %in% ereMAPs ~ T,
#                                         !(locus %in% ereMAPs) ~ F))

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, ereMAP == T & significant == T), aes(x = log2FoldChange, y = -log10(padj), fill = significant), size = 2, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, ereMAP == F), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#fb9a99')) +
  labs(fill= "") +
  ggrepel::geom_label_repel(data = filter(input, ereMAP == T & significant == T),
            aes(x = log2FoldChange,
                y = -log10(padj), 
                label = gene),
            force = 2)


  ggrepel::geom_text_repel(data = filter(input, ereMAP == T), 
                           aes(x = log2FoldChange, 
                               y = -log10(padj), 
                               label = locus))

## Overlap status annotated
  
input = as.data.frame(test)

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, significant == T), aes(x = log2FoldChange, y = -log10(padj), fill = overlap_expression), size = 2, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, significant == F), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#fb9a99')) +
  labs(fill= "")

## LIONS annotated

input = results_df_local_ERE
input = mutate(input, status = case_when(locus %in% LIONS_HI$TEid ~ 'LIONS_HI',
                                         locus %in% LIONS_LO$TEid ~ 'LIONS_LO',
                                         ID %in% ereMAP_loci ~ 'ereMAP',
                                        T ~ 'none'))

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, status != 'none' & significant == T), aes(x = log2FoldChange, y = -log10(padj), fill = status), size = 3, alpha = 1, shape = 21, stroke = 0) +
  geom_point(data = subset(input, status == 'none'), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  scale_fill_brewer(palette = 'Set1') +
  labs(fill= "")

#################################################################
# SalmonTE
#################################################################

input = Salmon_results_transcripts_df
input$TE = row.names(input)

volcano_plot = ggplot(data = input, aes(x = log2FoldChange, y = -log10(padj), fill = significant)) +
  geom_point(size = 2.5, alpha = 1, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  scale_fill_manual(values = c('#9B9A99', '#e41a1c')) +
  labs(fill= "") +
  ggrepel::geom_text_repel(data = filter(input, significant == T), aes(label = TE))

#################################################################
# TE_local + edgeR
#################################################################

input = edgeR_results_annotated

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = -logFC, y = -log10(FDR)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, ereMAP == T & significant == T), aes(x = -logFC, y = -log10(FDR), fill = significant), size = 2, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, ereMAP == F), aes(x = -logFC, y = -log10(FDR)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(FDR)')) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#fb9a99')) +
  labs(fill= "") +
  ggrepel::geom_label_repel(data = filter(input, ereMAP == T & significant == T),
                            aes(x = -logFC,
                                y = -log10(FDR), 
                                label = ID),
                            force = 2)

# Overlap

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = -logFC, y = -log10(FDR)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input,significant == T), aes(x = -logFC, y = -log10(FDR), fill = overlap_expression), size = 2, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(input, significant == F), aes(x = -logFC, y = -log10(FDR)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(FDR)')) +
  scale_fill_manual(values = c('#e41a1c', '#984ea3', '#4daf4a', '#fb9a99')) +
  labs(fill= "")

#################################################################
# Plot
#################################################################

volcano_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank(),
                                  legend.text = element_text(size = 15),
                                  legend.title = element_text(size = 18),
                                  legend.position = c(0.2, 0.93),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/volcano_TElocal.png", 
       width = 5.25, height = 5.25, units = "in")
