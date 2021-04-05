library(tidyverse)
library(ggplot2)
library(RColorBrewer)

## Calculate %TEs

df_gene_counts = data.frame(gene_counts = colSums(counts(dds_transcripts_gene, normalized = T)))
df_gene_counts = cbind(ID = rownames(df_gene_counts), df_gene_counts)

df_TE_counts = data.frame(TE_counts = colSums(counts(dds_transcripts_TE, normalized = T)))
df_TE_counts = cbind(ID = rownames(df_TE_counts), df_TE_counts)

df = merge(df_gene_counts, df_TE_counts, by = 'ID')

df = mutate(df, total = gene_counts + TE_counts)
df = mutate(df, percent_TEs = (TE_counts / total)*100)

df = separate(df, col = ID, into = c('patient', 'tissue'), sep = '_')

df = mutate(df, tissue = fct_reorder(tissue, percent_TEs))

## Boxplot

plot = ggplot(data = df, aes(x = tissue, y = percent_TEs)) +
  geom_boxplot(fill = "#9ECAE1", width = 0.6) +
  geom_point() +
  ylab('% of reads mapping to TEs') +
  xlab('Tissue/cell')

plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 16),
                            plot.subtitle = element_text(size = 12),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_line(color = 'gray'),
                            axis.text.x = element_text(size = 13, vjust = -.5),
                            axis.text.y = element_text(size = 13),
                            axis.title.x = element_text(size = 15, margin = margin(t = 12.5)),
                            axis.title.y = element_text(size = 15, margin = margin(r = 7.5)),
                            axis.line = element_line(size = 0.8),
                            panel.border = element_blank(),
                            legend.text = element_text(size = 12),
                            legend.title = element_text(size = 14))

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/overall_TE_abundance.png", 
       width = 7, height = 4, units = "in")

## One-way ANOVA

aov_res = aov(percent_TEs ~ tissue, data = df)
summary(aov_res)
tukey_res = TukeyHSD(aov_res)
