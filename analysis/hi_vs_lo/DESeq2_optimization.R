library(ggplot2)
library(dplyr)

gene = readRDS('~/TE_thymus/analysis/cluster/objects/sig_diff_gene_count.rds')
TE = readRDS('~/TE_thymus/analysis/cluster/objects/sig_diff_TE_count.rds')

gene = as.data.frame(gene)
TE = as.data.frame(TE)

gene = cbind(independent_filtering = rownames(gene), gene)
TE = cbind(independent_filtering = rownames(TE), TE)

gene = pivot_longer(data = gene, cols = '2':'20', names_to = 'min_reads', values_to = 'gene_sigdiff_count')
TE = pivot_longer(data = TE, cols = '2':'20', names_to = 'min_reads', values_to = 'TE_sigdiff_count')

sigdiff_count = merge(gene, TE)

sigdiff_count$min_reads = factor(sigdiff_count$min_reads, levels = c('2', '5', '10', '20'))

## Gene expression

plot = ggplot(data = sigdiff_count, aes(x = min_reads, y = gene_sigdiff_count, fill = independent_filtering)) +
  geom_bar(color = 'black', stat = 'identity', position = 'dodge') +
  scale_fill_brewer(palette = "Set1") +
  xlab('Minimum normalized reads per gene') +
  ylab('Number of differentially expressed genes') +
  ggtitle('Effect of DESeq2 parameters on detection of differential expression', 'Gene expression')

## TE expression

plot = ggplot(data = sigdiff_count, aes(x = min_reads, y = TE_sigdiff_count, fill = independent_filtering)) +
  geom_bar(color = 'black', stat = 'identity', position = 'dodge') +
  scale_fill_brewer(palette = "Set1") +
  xlab('Minimum normalized reads per copy') +
  ylab('Number of differentially expressed copies') +
  ggtitle('Effect of DESeq2 parameters on detection of differential expression', 'TE expression') +
  scale_y_continuous(limits = c(0,15000), breaks = c(0, 2500, 5000, 7500, 10000, 12500, 15000))

plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 16),
                   plot.subtitle = element_text(size = 14),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_line(color = 'grey'),
                   panel.grid.minor.y = element_blank(),
                   axis.text.x = element_text(size = 13),
                   axis.text.y = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   axis.line = element_line(size = 0.8),
                   panel.border = element_blank(),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 14))

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/21-03-25/DESeq2_optimization_gene_count.png", 
       width = 24, height = 12, units = "cm")