library(ggplot2)
library(dplyr)
library(tidyr)

## Preparing the input data


input = results_df_local_TE
input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))
input = filter(input, significant == TRUE)

input = mutate(input, up_regulated = case_when(log2FoldChange > 0 ~ T,
                                               log2FoldChange < 0 ~ F))

input = mutate(input, detected_by_TE = case_when(gene %in% results_df_transcripts_TE_sigdiff$gene ~ T,
                                                 !(gene %in% results_df_transcripts_TE_sigdiff$gene) ~ F))

#################################################################
# Plotting changes in expression at the element level
#################################################################

## Take a random subset

input_random_subset = sample_n(input, 200)

## Plot fold change

plot = ggplot(data = input_random_subset, aes(x = gene, y = log2FoldChange, color = detected_by_TE)) +
  geom_bar(stat = 'summary', fun.y = 'mean') + 
  geom_point() +
  xlab('Element') +
  ylab(expression('Log'[2]*' Fold Change')) +
  ggtitle('TE expression in mTEC-HI vs mTEC-LO', 'Differentially regulated copies only (random subset of elements)')


plot + theme(axis.text.x = element_text(angle = 90, size = 10),
             plot.title = element_text(face = 'bold', size = 20),
             plot.subtitle = element_text(size = 14))

## Plot signed p value

input_random_subset = mutate(input_random_subset, signed_log10padj = -log10(padj) * sign(log2FoldChange))

ggplot(data = input, aes(x = gene, y = signed_log10padj)) +
  geom_point(alpha = 0.6, size = 0.5) +
  facet_grid(. ~ class)




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


