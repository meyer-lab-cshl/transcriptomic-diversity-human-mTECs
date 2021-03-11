#################################################################
# Volcano
#################################################################

input = filter(df, Class == 'protein_coding')

## Coloured by significance (statistical and biological)

gene_list = 'AIRE'

volcano_plot = ggplot(data = input, aes(x = log2FoldChange, y = -log10(padj), colour = overall_significant)) +
  geom_point(alpha = 0.6, aes(colour = overall_significant)) +
  xlab(expression('Log'[2]*' FC')) +
  ylab(expression('-Log'[10]*' P value')) +
  ggtitle('mTEC-hi vs mTEC-lo', 'Protein-coding gene expression') +
  scale_colour_manual(values = c('#9B9A99', "red")) +
  geom_label_repel(
    data = subset(input, GeneSymbol == gene_list),
    aes(label = subset(input, GeneSymbol == gene_list)$GeneSymbol),
    size = 3.5, 
    box.padding = unit(0.9, 'lines'), 
    point.padding = unit(0.2, 'lines'),
    max.overlaps = 20)

#xlim(-3, 3) +
#geom_vline(xintercept = c(1, -1), linetype = 'dashed') +
#geom_hline(yintercept = -log10(0.05), linetype = 'dashed')

## Plot

volcano_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_gene_volcano_plot.png", 
       width = 20, height = 15, units = "cm")