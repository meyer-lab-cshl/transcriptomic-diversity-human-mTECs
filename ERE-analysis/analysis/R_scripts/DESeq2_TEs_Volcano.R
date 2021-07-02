library(ggplot2)
library(ggrepel)
library(dplyr)

#################################################################
# TE local
#################################################################

input = results_df_local_TE

input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))

## Coloured by significance (statistical and biological)

volcano_plot = ggplot(data = input, aes(x = log2FoldChange, y = -log10(padj), colour = significant)) +
  geom_point(alpha = 0.6, aes(colour = significant)) +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  scale_colour_manual(values = c('#9B9A99', "red")) +
  geom_hline(yintercept = -log10(0.1), linetype = 'dashed') + 
  labs(color= "FDR < 0.1")

#geom_vline(xintercept = c(1, -1), linetype = 'dashed') +
#geom_hline(yintercept = -log10(0.05), linetype = 'dashed')

## Colored by TE class

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = class), size = 1, alpha = 0.8, shape = 21, stroke = 0.1) +
  geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value'))

#  scale_colour_manual(values = c('#e41a1c', "#377eb8", "#4daf4a", "#984ea3"))

## Colored by differentially regulated in TE_transcripts

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = diff_reg_both, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('red', 0.6)) +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  ggtitle('mTEC-hi vs mTEC-lo', 'TE expression: TElocal')

## Coloured by spatial overlap with DE genes 

volcano_plot = ggplot(data = subset(input, overlap_status != '<NA>'), aes(x = log2FoldChange, y = -log10(padj), colour = overlap_status)) +
  geom_point(alpha = 0.6) +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  ggtitle('mTEC-hi vs mTEC-lo', 'TE expression: TE local') +
  ylim(0, 90)
  scale_colour_manual(values = c('#9B9A99', "red"))

## Plot

volcano_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank(),
                                  legend.text = element_text(size = 11),
                                  legend.title = element_text(size = 14))

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/local_volcano.png", 
       width = 11, height = 5.25, units = "in")

#################################################################
# TE transcripts
#################################################################

input = results_df_transcripts_TE

input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))

## Coloured by significance (statistical and biological)

volcano_plot = ggplot(data = input, aes(x = log2FoldChange, y = -log10(padj), colour = significant)) +
  geom_point(alpha = 0.6, aes(colour = significant)) +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  ggtitle('mTEC-hi vs mTEC-lo', 'TE expression') +
  scale_colour_manual(values = c('#9B9A99', "red")) + 
  geom_hline(yintercept = -log10(0.1), linetype = 'dashed')




  guides(colour = FALSE) +
  geom_label_repel(
    data = subset(input, overall_significant == TRUE),
    aes(label = subset(input, overall_significant == TRUE)$family),
    size = 3.5, 
    box.padding = unit(0.9, 'lines'), 
    point.padding = unit(0.2, 'lines'),
    max.overlaps = 20)

#geom_vline(xintercept = c(1, -1), linetype = 'dashed') +
#geom_hline(yintercept = -log10(0.05), linetype = 'dashed')

## Colored by TE class

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = class), size = 2.5, alpha = 0.8, shape = 21, stroke = 0.4) +
  geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  xlim(-2, 3) +
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
                                  legend.title = element_text(size = 14))

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/transcripts_volcano.png", 
       width = 11, height = 5.25, units = "in")
