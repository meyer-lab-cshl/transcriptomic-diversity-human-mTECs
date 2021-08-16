library(DESeq2)
library(tidyr)
library(ggplot2)
library(glue)

functions_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_functions/"
functions = c('extract_subset', 'differential_expression', 'process_DESeq2_results', 'build_count_table', 'generate_contingency')
f
for (i in functions){
  
  load(glue('{functions_directory}{i}'))
  
}

#################################################################
# Stacked bars
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

#################################################################
#  (supplement B?)
#################################################################

count_table = build_count_table(dds_transcripts_TE, 
                                results_df_transcripts_TE, 
                                group = c('all', 'down_regulated', 'up_regulated'),
                                mode = 'LTR_family',
                                by = 'normalized_reads')

count_table = mutate(count_table, grouped_class = case_when(family == 'ERV1' ~ 'ERV1',
                                                            family == 'ERVK' ~ 'ERVK',
                                                            family == 'ERVL' ~ 'ERVL',
                                                            family == 'ERVL-MaLR' ~ 'ERVL-MaLR',
                                                            family == 'Gypsy' ~ 'Other',
                                                            family == 'LTR' ~ 'Other'))

count_table = group_by(count_table, group, grouped_class) %>% 
  summarize(percent = sum(percent))

count_table$grouped_class = factor(count_table$grouped_class, levels = c('Other', 'ERVL-MaLR', 'ERVL', 'ERVK', 'ERV1'))

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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/SB_stacked-bars.png", 
       width = 6, height = 5, units = "in")


#################################################################
# Odds ratio
#################################################################

input = results_df_transcripts_ERE 

class = vector()
p_value = vector()
odds_ratio = vector()
lower_interval = vector()
upper_interval = vector()

for (i in 1:length(unique(input$class))){
  
  query_class = unique(input$class)[i]
  
  column_1 = c(nrow(subset(input, significant == T & log2FoldChange > 0 & class == query_class)), 
               nrow(subset(input, !(significant == T & log2FoldChange > 0) & class == query_class)))
  
  column_2 = c(nrow(subset(input, significant == T & log2FoldChange > 0 & class != query_class)), 
               nrow(subset(input, !(significant == T & log2FoldChange > 0) & class != query_class)))
  
  contingency = data.frame(class = column_1, 'Not' = column_2, row.names = c('Up', 'Not up'))
  
  print(query_class)
  print(contingency)
  
  class[i] = query_class
  p_value[i] = fisher.test(x = contingency)[[1]]
  odds_ratio[i] = fisher.test(x = contingency)[[3]]
  lower_interval[i] = fisher.test(x = contingency)$conf.int[1]
  upper_interval[i] = fisher.test(x = contingency)$conf.int[2]
  
}

p_value = p.adjust(p_value, method = 'bonferroni')

output = data.frame(class = class, 
                    p_value = p_value, 
                    odds_ratio = odds_ratio, 
                    lower_interval = lower_interval, 
                    upper_interval = upper_interval) %>%
  mutate(significant = case_when(p_value < 0.001 ~ '***', 
                                 p_value < 0.01 ~ '**',
                                 p_value < 0.05 ~ '*',
                                 p_value >= 0.05 ~ '')) %>%
  mutate(class = forcats::fct_reorder(class, odds_ratio))

odds_ratio_plot = ggplot(data = output, aes(x = class, y = odds_ratio)) + 
                  geom_point() +
                  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_errorbar(aes(ymin = lower_interval, ymax = upper_interval), width = 0) +
  xlab('') +
  ylab('Odds ratio') +
  geom_text(aes(label = significant), nudge_y = 1, size = 6)
  

odds_ratio_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_text(size = 15, margin = margin(r = 7.5)),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank(),
                                  legend.text = element_text(size = 15),
                                  legend.title = element_text(size = 18),
                                  legend.position = c(0.2, 0.93),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/odds_ratio.png", 
       width = 5, height = 5, units = "in")
  
#################################################################
# TREs
#################################################################

TREs = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/TREs_Larouche.csv', header = T) %>%
  select(c(ERE.family, Tissues.overexpressing))

results_df_transcripts_ERE = mutate(results_df_transcripts_ERE, TRE = case_when(gene %in% TREs$ERE.family ~ T,
                                                                                !(gene %in% TREs$ERE.family) ~ F))

generate_contingency = function(input, condition_A, condition_B){
  
  col_1 = vector()
  col_2 = vector()
  frequency = vector()
  index = 1
  
  for (a in 1:length(condition_A)){
    
    for (b in 1:length(condition_B)){
      
      col_1[index] = names(condition_A)[a]
      col_2[index] = names(condition_B)[b]
      
      frequency[index] = nrow(subset(input, eval(parse(text=condition_A[a])) & eval(parse(text=condition_B[b]))))
      
      index = index + 1
      
    }
    
  }
  
  count_table  = data.frame(condition_A = col_1,
                      condition_B = col_2,
                      frequency = frequency)
  
  output = xtabs(frequency ~ condition_A+condition_B, data=count_table)
  
  return(output)
  
}

output = generate_contingency(input = results_df_transcripts_ERE,
                     condition_A = list('up_TEs' = 'significant == T & log2FoldChange > 0',
                                       'down_TEs' = 'significant == T & log2FoldChange < 0',
                                       'unchanged_TEs' = 'significant == F'),
                     condition_B = list('TRE' = 'TRE == T', 'not_TRE' = 'TRE == F'))

vcd::mosaic(~condition_A+condition_B, data = output, direction = c('v', 'h'), shade = T)

#################################################################
# Class
#################################################################

output = generate_contingency(input = results_df_transcripts_ERE,
                              condition_A = list('up_TEs' = 'significant == T & log2FoldChange > 0',
                                                 'down_TEs' = 'significant == T & log2FoldChange < 0',
                                                 'unchanged_TEs' = 'significant == F'),
                              condition_B = list('LTR' = 'class == \'LTR\'', 
                                                 'SINE' = 'class == \'SINE\'',
                                                 'LINE' = 'class == \'LINE\''))

names(dimnames(output)) <- c("Expression_status", "ERE_class")


png("~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/mosaic_expression_vs_class.png", width = 400, height = 400)
vcd::mosaic(~Expression_status+ERE_class, data = output, 
            direction = c('v', 'h'), 
            shade = T,
            set_labels = list(Expression_status = c('down_TEs' = 'Down',
                                                    'unchanged_TEs' = 'Unchanged',
                                                    'up_TEs' = 'Up')))
dev.off()

#################################################################
# LIONS
#################################################################

lions = read.csv('~/Desktop/thymus-epitope-mapping/ERE-analysis/LIONS/pt214_hi_1.lion', header = T, sep = '\t')

TE_initiations = tidyr::separate(data = lions, col = repeatName, sep = ':', into = c('gene', 'class', 'family'))$gene

subset(results_df_transcripts_ERE, gene %in% TE_initiations & significant == T)
  