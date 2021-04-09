library(ggplot2)
library(dplyr)
library(tidyr)

## Preparing the input data


input = results_df_local_TE
input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))
#input = filter(input, significant == TRUE)

#input = mutate(input, up_regulated = case_when(log2FoldChange > 0 ~ T,
#                                               log2FoldChange < 0 ~ F))

input = mutate(input, detected_by_TE = case_when(gene %in% results_df_transcripts_TE_sigdiff$gene ~ T,
                                                 !(gene %in% results_df_transcripts_TE_sigdiff$gene) ~ F))

input = mutate(input, log2FoldChange_transcripts = )

#################################################################
# Plotting changes in expression at the element level
#################################################################

## Take a random subset

random_element_subset = sample(results_df_transcripts_TE$gene, 30)

local_input = filter(results_df_local_TE, gene %in% random_element_subset)

transcripts_input = filter(results_df_transcripts_TE, gene %in% random_element_subset) %>% 
  arrange(log2FoldChange) %>%
  mutate(gene = factor(gene, levels = gene))

plot = ggplot(data = transcripts_input, aes(x = gene, y = log2FoldChange, fill = significant)) +
  geom_bar(stat = 'identity') + 
  geom_jitter(data = local_input, aes(x = gene, y = log2FoldChange, color = significant), size = 0.5, alpha = 0.6, width = 0.3) +
  xlab('Element') +
  ylab(expression('Log'[2]*' Fold Change')) 

  
##ggtitle('TE expression in mTEC-HI vs mTEC-LO', 'TEtranscripts and TElocal compared (random subset of elements)')

plot = ggplot(data = transcripts_input, aes(x = gene, y = log2FoldChange)) +
  geom_bar(stat = 'identity', aes(fill = significant)) + 
  geom_violin(data = local_input, aes(x = gene, y = log2FoldChange)) +
  xlab('Element') +
  ylab(expression('Log'[2]*' Fold Change')) +
  ggtitle('TE expression in mTEC-HI vs mTEC-LO', 'TEtranscripts and TElocal compared (random subset of elements)')

## Only elements detected by TE_transcirpts

local_input = filter(results_df_local_TE, gene %in% results_df_transcripts_TE_sigdiff$gene)

transcripts_input = filter(results_df_transcripts_TE, gene %in% results_df_transcripts_TE_sigdiff$gene) %>% 
  arrange(log2FoldChange) %>%
  mutate(gene = factor(gene, levels = gene))

plot = ggplot(data = transcripts_input, aes(x = gene, y = log2FoldChange, fill = significant)) +
  geom_bar(stat = 'identity') + 
  geom_point(data = local_input, aes(x = gene, y = log2FoldChange, color = significant), size = 0.5, alpha = 0.6) +
  xlab('Element') +
  ylab(expression('Log'[2]*' Fold Change')) +
  ggtitle('TE expression in mTEC-HI vs mTEC-LO', 'TEtranscripts and TElocal compared (elements detected by TEtranscripts only)')

plot = ggplot(data = transcripts_input, aes(x = gene, y = log2FoldChange)) +
  geom_bar(stat = 'identity', aes(fill = significant)) + 
  geom_violin(data = local_input, aes(x = gene, y = log2FoldChange)) +
  xlab('Element') +
  ylab(expression('Log'[2]*' Fold Change')) +
  ggtitle('TE expression in mTEC-HI vs mTEC-LO', 'TEtranscripts and TElocal compared (random subset of elements)')

## Plot

plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 16),
                          plot.subtitle = element_text(size = 12),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          panel.grid.major.y = element_line(color = 'gray'),
                          panel.grid.minor.y = element_blank(),
                          axis.text.x = element_text(size = 12, angle = 90),
                          axis.text.y = element_text(size = 13),
                          axis.title.x = element_text(size = 15, margin = margin(t = 12.5)),
                          axis.title.y = element_text(size = 15, margin = margin(r = 7.5)),
                          axis.line = element_line(size = 0.8),
                          panel.border = element_blank(),
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size = 14),
                          axis.ticks.length = unit(.15, 'cm'))

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/local_vs_transcripts_violin.png", 
       width = 11.5, height = 6, units = "in")

#################################################################
# Only elements detected by TE_transcripts
#################################################################




#################################################################
# 
#################################################################

output = input %>% group_by(gene, up_regulated, detected_by_TE) %>% summarize(total = n())
output = spread(data = output, key = up_regulated, value = total)
output$down = replace_na(output$'FALSE', 0)
output$up = replace_na(output$'TRUE', 0)
output = select(output, -c('FALSE', 'TRUE'))
output = mutate(output, total = up + down)
output = mutate(output, percent_up = (up / total) * 100)

output = mutate(output, consistent = case_when((percent_up == 100) | (percent_up == 0) ~ TRUE,
                                               (percent_up != 100) & (percent_up != 0) ~ FALSE))

output = output %>% group_by(detected_by_TE, consistent) %>% summarize(count = n())

output = as.data.frame(output)

ggplot(data = output, aes(x = consistent, y = count, fill = detected_by_TE)) +
  geom_col(colour = 'black', position = 'fill')

## Fisher's exact test

contingency_table = matrix(c(65, 395, 18, 179), 
                          nrow = 2,
                          dimnames = list('detected_by_transcripts' = c('True', 'False'), 'consistent' = c('True', 'False')))

fisher.test(contingency_table, alternative = 'greater')




###

output = as.data.frame(output)

ggplot(data = output, aes(x = reorder(gene, -percent_up), y = percent_up, fill = detected_by_TE)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = "Set1")


