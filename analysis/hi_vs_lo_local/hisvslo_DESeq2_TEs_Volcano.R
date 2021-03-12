library(ggplot2)
library(ggrepel)

#################################################################
# Volcano
#################################################################

input = results_df

input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))

## Coloured by significance (statistical and biological)

volcano_plot = ggplot(data = input, aes(x = log2FoldChange, y = -log10(padj), colour = overall_significant)) +
  geom_point(alpha = 0.6, aes(colour = overall_significant)) +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  ggtitle('mTEC-hi vs mTEC-lo', 'TE expression') +
  scale_colour_manual(values = c('#9B9A99', "red")) +
  guides(colour = FALSE) +
  geom_label_repel(
    data = subset(results_df, overall_significant == TRUE),
    aes(label = subset(results_df, overall_significant == TRUE)$class),
    size = 3.5, 
    box.padding = unit(0.9, 'lines'), 
    point.padding = unit(0.2, 'lines'),
    max.overlaps = 20)

#geom_vline(xintercept = c(1, -1), linetype = 'dashed') +
#geom_hline(yintercept = -log10(0.05), linetype = 'dashed')

## Colored by TE class

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, overall_significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), color = class), alpha = 0.6) +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  ggtitle('mTEC-hi vs mTEC-lo', 'TE expression')

#  scale_colour_manual(values = c('#e41a1c', "#377eb8", "#4daf4a", "#984ea3"))

## Plot

volcano_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo_local/Plots/hi_vs_lo_TEs_volcano_plot.png", 
       width = 20, height = 15, units = "cm")
