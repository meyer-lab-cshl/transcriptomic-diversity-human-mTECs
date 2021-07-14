library(DESeq2)
library(tidyr)
library(ggplot2)
library(glue)

functions_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_functions/"
functions = c('extract_subset', 'differential_expression', 'process_DESeq2_results', 'build_count_table')

for (i in functions){
  
  load(glue('{functions_directory}{i}'))
  
}

#################################################################
# DESeq2
#################################################################

count_table_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/"
data = read.table(glue('{count_table_directory}TE_transcripts_hi_vs_lo.cntTable'),header=T,row.names=1)

## Run DESeq2

TE_data = extract_subset(mode = 'TE', input = data)
dds_transcripts_TE = differential_expression(TE_data, design=~patient + tissue)

vs_dds_transcripts_TE = vst(dds_transcripts_TE, blind=FALSE)

results_transcripts_TE = results(dds_transcripts_TE, 
                                 contrast = c('tissue', 'mTEC-hi', 'mTEC-lo'), 
                                 independentFiltering = F)
results_df_transcripts_TE = process_DESeq2_results(results = results_transcripts_TE, mode = 'TE_transcripts')
results_df_transcripts_TE_sigdiff = filter(results_df_transcripts_TE, significant == T)

#################################################################
# Volcano plot (B)
#################################################################

## Prepare the input

input = results_df_transcripts_TE

input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))

input = mutate(input, grouped_class = case_when(class == 'LINE' ~ 'LINE',
                                                   class == 'LTR' ~ 'LTR',
                                                   class == 'SINE' ~ 'SINE',
                                                   class == 'Retroposon' ~ 'Other',
                                                   class == 'Satellite' ~ 'Other',
                                                   class == 'RC' ~ 'Other',
                                                   class == 'DNA' ~ 'DNA',
                                                   class == 'RNA' ~ 'Other',
                                                   class == 'Unknown' ~ 'Other'))

input$grouped_class = factor(input$grouped_class, levels = c('LTR', 'DNA', 'LINE', 'SINE', 'Other'))

## Plot

volcano_plot = ggplot() +
  geom_point(data = input, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = grouped_class), size = 2.5, alpha = 1, shape = 21, stroke = 0) +
  geom_point(data = subset(input, significant == FALSE), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 1, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  xlim(-2, 3) +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff', '#55A257', '#E93C00', '#9A9A9A')) +
  labs(fill= "")

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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/B_volcano_v2.png", 
       width = 5.25, height = 5.25, units = "in")

#################################################################
# Stacked bars (C)
#################################################################

count_table = build_count_table(dds_transcripts_TE, 
                                results_df_transcripts_TE, 
                                group = c('all', 'down_regulated', 'up_regulated'),
                                mode = 'class',
                                by = 'normalized_reads')

count_table = mutate(count_table, grouped_class = case_when(class == 'LINE' ~ 'LINE',
                                                            class == 'LTR' ~ 'LTR',
                                                            class == 'SINE' ~ 'SINE',
                                                            class == 'Retroposon' ~ 'Other',
                                                            class == 'Satellite' ~ 'Other',
                                                            class == 'RC' ~ 'Other',
                                                            class == 'DNA' ~ 'DNA',
                                                            class == 'RNA' ~ 'Other',
                                                            class == 'Unknown' ~ 'Other'))


count_table = group_by(count_table, group, grouped_class) %>% 
  summarize(percent = sum(percent))

count_table$grouped_class = factor(count_table$grouped_class, levels = c('Other', 'SINE', 'LINE', 'DNA', 'LTR'))

bar_chart = ggplot(count_table, aes(x = group, y = percent, fill = grouped_class)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  xlab('') +
  ylab('Fraction of normalized reads') +
  labs(fill= "") +
  scale_x_discrete(labels = c('All', 'Downregulated', 'Upregulated')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c('#9A9A9A', '#E93C00', '#55A257', '#dd8452ff', '#4c72b0ff'))

bar_chart + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                               plot.subtitle = element_text(size = 14),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text.x = element_text(size = 13, margin = margin(t = 6)),
                               axis.text.y = element_text(size = 14),
                               axis.title.y = element_text(size = 14),
                               axis.title.x = element_text(size = 14, margin = margin(t = 6)),
                               axis.line = element_line(size = 0.8),
                               panel.border = element_blank(),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 14),
                               legend.position="top")

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/C_stacked-bars.png", 
       width = 5, height = 5, units = "in")

## Chi-square test

contingency = pivot_wider(count_table, id_cols = grouped_class, names_from = group, values_from = sum) %>% 
  mutate_if(is.numeric, list(~replace_na(., 0))) %>% 
  as.data.frame()

rownames(contingency) = contingency$grouped_class
contingency = dplyr::select(contingency, -grouped_class)

chisq.test(contingency)

chisq.posthoc.test::chisq.posthoc.test(contingency)
