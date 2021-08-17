library(tidyverse)
library(GenomicRanges)
library(glue)

working_directory = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis'

## Import required variables
GRanges_gene_extended = readRDS(file = glue('{working_directory}/R_variables/GRanges_gene_extended'))
GRanges_ERE_start = readRDS(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_ERE_start')
up_genes = readRDS(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/up_genes')
down_genes = readRDS(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/down_genes')

##Import required functions

functions = c('generate_contingency')

for (i in functions){
  
  load(glue('{working_directory}/R_functions/{i}'))
  
}

#################################################################
# Overlap analysis
#################################################################

find_frequency_of_overlaps = function(gene_universe, TE_universe, gene_sets, TE_sets){
  
  for (a in 1:length(TE_sets)){
    
    total_overlaps = length(findOverlaps(query = TE_sets[[a]],
                                         subject = gene_universe))
    
    print(total_overlaps)
    
    percent = vector()
    for (b in 1:length(gene_sets)){
      
      subset_overlap = length(findOverlaps(query = TE_sets[[a]],
                                           subject = gene_sets[[b]]))
      
      percent[b] = subset_overlap / total_overlaps * 100
      
    }
    
    table_entry = data.frame(gene_set = names(gene_sets),
                             percent_overlap = percent,
                             TE_set = names(TE_sets)[a])
    
    if (a == 1){
      
      frequency_table = table_entry
      
    }
    
    else{
      
      frequency_table = bind_rows(frequency_table, table_entry)
      
    }
    
  }
  
  ## Reformat the frequency table to a contingency table
  
  contingency_table = xtabs(percent_overlap ~ gene_set+TE_set, data=frequency_table)
  
  output = list('frequency_table' = frequency_table,
                'contingency_table' = contingency_table)
  
  return(output)
  
}

## Differentially regulated genes

TE_sets = list('up_EREs' = subset(GRanges_ERE_start, locus %in% up_EREs),
               'unchanged_EREs' = subset(GRanges_ERE_start, !(locus %in% up_EREs) & !(locus %in% down_EREs)),
               'down_EREs' = subset(GRanges_ERE_start, locus %in% down_EREs))

gene_sets = list('up_genes' = subset(GRanges_gene_extended, Geneid %in% up_genes),
               'unchanged_genes' = subset(GRanges_gene_extended, !(Geneid %in% up_genes) & !(Geneid %in% down_genes)),
               'down_genes' = subset(GRanges_gene_extended, Geneid %in% down_genes))

final_output = find_frequency_of_overlaps(gene_universe = GRanges_gene_extended,
                                          TE_universe = GRanges_ERE_start,
                                          gene_sets = gene_sets,
                                          TE_sets = TE_sets)

vcd::mosaic(~TE_set+gene_set, data = final_output$contingency_table, 
            direction = c('v', 'h'), 
            shade = T,
            labeling = vcd::labeling_border(rot_labels = c(0,0,0,0),
                                            offset_varnames = c(1,0,0,7),
                                            just_labels=c('center', 'right')),
            margins = c(0,1,0,5))

## Genes of interest

gene_sets = list('AIRE_genes' = subset(GRanges_gene_extended, Geneid %in% AIRE_genes),
                 'FEZF2_genes' = subset(GRanges_gene_extended, Geneid %in% FEZF2_genes),
                 'Housekeeping_genes' = subset(GRanges_gene_extended, Geneid %in% housekeeping_genes), 
                 'other_genes' = subset(GRanges_gene_extended, !(Geneid %in% AIRE_genes) & !(Geneid %in% FEZF2_genes) & !(Geneid %in% housekeeping_genes)))

gene_sets = list('lncRNA' = subset(GRanges_gene_extended, Class == 'lncRNA'),
                 'protein' = subset(GRanges_gene_extended, Class == 'protein_coding'),
                 'other' = subset(GRanges_gene_extended, Class != 'lncRNA' & Class != 'protein_coding'))

final_output = find_frequency_of_overlaps(gene_universe = GRanges_gene_extended,
                                          TE_universe = GRanges_ERE_start,
                                          gene_sets = gene_sets,
                                          TE_sets = TE_sets)

vcd::mosaic(~TE_set+gene_set, data = final_output$contingency_table, 
            direction = c('v', 'h'), 
            shade = T,
            labeling = vcd::labeling_border(rot_labels = c(0,0,0,0),
                                            offset_varnames = c(1,0,0,8),
                                            just_labels=c('center', 'right')),
            margins = c(0,1,0,5))



# Plot bar chart

r = '#e41a1c'
p = '#984ea3'
g = '#4daf4a'

input = final_output$frequency_table
input$class = factor(input$class, levels = c('LTR', 'SINE', 'LINE', 'Retroposon'))

bar_chart = ggplot(input, aes(x = TE_set, y = percent_overlap, fill = gene_set)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  xlab('') +
  ylab('Fraction of overlap events') +
  labs(fill= "") +
  scale_x_discrete(labels = c('Down EREs', 'Unchanged EREs', 'Up EREs')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  guides(fill = guide_legend(reverse = T)) +
  scale_fill_manual(labels = c('Down genes', 'Unchanged genes', 'Up genes'), values = c(p, g, r))

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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/fraction_of_overlaps_with_GOIs.png", 
       width = 5, height = 5, units = "in")

################################################################################
# Overlap between TEs and expressed genes
################################################################################

results_df_local_ERE = readRDS(file = glue::glue('{working_directory}/R_variables/results_df_local_ERE'))
#overlap_annotated_GRanges_gene_extended = readRDS(file = glue::glue('{working_directory}/R_variables/overlap_annotated_GRanges_gene_extended'))
overlap_annotated_GRanges_gene_extended = input

input = merge(subset(results_df_local_ERE, significant == T), as.data.frame(overlap_annotated_GRanges_gene_extended), by = 'locus')%>%
  mutate(TE_expression = case_when(significant == T & log2FoldChange > 0 ~ 'up',
                                   significant == T & log2FoldChange < 0 ~ 'down',
                                   T ~ 'unchanged'))

condition_A = list('unchanged' = 'overlap_expression == "unchanged"',
                   'up' = 'overlap_expression == "up"',
                   'down' = 'overlap_expression == "down"',
                   'none' = 'overlap_expression == "none"')

condition_B = list('TE_up' = 'TE_expression == "up"',
                   'TE_down' = 'TE_expression == "down"')

output = generate_contingency(input = input, 
                     condition_A = condition_A, 
                     condition_B = condition_B)

vcd::mosaic(~condition_B+condition_A, data = output, 
            direction = c('v', 'h'), 
            shade = T)

## Fishers exact

col_1 = c(nrow(subset(input, TE_expression == 'up' & overlap_expression == 'down')), 
          nrow(subset(input, TE_expression == 'up' & overlap_expression != 'down')))

col_2 = c(nrow(subset(input, TE_expression == 'down' & overlap_expression == 'down')), 
          nrow(subset(input, TE_expression == 'down' & overlap_expression != 'down')))

contingency = data.frame(col_1 = col_1, col_2 = col_2)

fisher.test(x = contingency)

expression = vector()
p_value = vector()
odds_ratio = vector()
lower_interval = vector()
upper_interval = vector()

for (i in 1:length(unique(input$overlap_expression))){
  
  query_class = unique(input$overlap_expression)[i]
  
  column_1 = c(nrow(subset(input, TE_expression == 'up' & overlap_expression == query_class)), 
               nrow(subset(input, TE_expression == 'down' & overlap_expression == query_class)))
  
  column_2 = c(nrow(subset(input, TE_expression == 'up' & overlap_expression != query_class)), 
               nrow(subset(input, TE_expression == 'down' & overlap_expression != query_class)))
  
  contingency = data.frame(class = column_1, 'Not' = column_2, row.names = c('Up', 'Down'))
  
  print(query_class)
  print(contingency)
  
  expression[i] = query_class
  p_value[i] = fisher.test(x = contingency)[[1]]
  odds_ratio[i] = fisher.test(x = contingency)[[3]]
  lower_interval[i] = fisher.test(x = contingency)$conf.int[1]
  upper_interval[i] = fisher.test(x = contingency)$conf.int[2]
  
}

p_value = p.adjust(p_value, method = 'bonferroni')

output = data.frame(class = expression, 
                    p_value = p_value, 
                    odds_ratio = odds_ratio, 
                    lower_interval = lower_interval, 
                    upper_interval = upper_interval) %>%
  mutate(significant = case_when(p_value < 0.001 ~ '***', 
                                 p_value < 0.01 ~ '**',
                                 p_value < 0.05 ~ '*',
                                 p_value >= 0.05 ~ '')) %>%
  mutate(class = forcats::fct_reorder(expression, odds_ratio))

odds_ratio_plot = ggplot(data = output, aes(x = expression, y = (odds_ratio))) + 
  geom_point() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_errorbar(aes(ymin = lower_interval, ymax = upper_interval), width = 0) +
  xlab('') +
  ylab('Odds ratio') +
  geom_text(aes(label = significant), nudge_y = 1, size = 6) +
  scale_y_continuous(trans='log2')


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