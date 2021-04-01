library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

df_gene_counts = data.frame(gene_counts = colSums(assay(dds_transcripts_gene)))
df_gene_counts = cbind(ID = rownames(df_gene_counts), df_gene_counts)

df_TE_counts = data.frame(TE_counts = colSums(assay(dds_transcripts_TE)))
df_TE_counts = cbind(ID = rownames(df_TE_counts), df_TE_counts)

df = merge(df_gene_counts, df_TE_counts, by = 'ID')

df = mutate(df, total = gene_counts + TE_counts)
df = mutate(df, percent_TEs = (TE_counts / total)*100)

df = separate(df, col = ID, into = c('patient', 'tissue'), sep = '_')

plot = ggplot(data = df, aes(x = tissue, y = percent_TEs)) +
  geom_boxplot(fill = "#9ECAE1") +
  geom_jitter(width = 0.1) +
  ylab('% of reads mapping to TEs') +
  xlab('Tissue/cell')

plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 16),
                            plot.subtitle = element_text(size = 12),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            axis.text.x = element_text(size = 12),
                            axis.text.y = element_text(size = 12),
                            axis.title = element_text(size = 12),
                            axis.line = element_line(size = 0.8),
                            panel.border = element_blank(),
                            legend.text = element_text(size = 12),
                            legend.title = element_text(size = 14))

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/overall_TE_abundance.png", 
       width = 11.61, height = 5.72, units = "in")
