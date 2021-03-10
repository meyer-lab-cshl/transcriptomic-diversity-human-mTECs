library(plyr)
library(dplyr)

#################################################################
# Enrichment 
#################################################################

calculate_enrichment_factor = function(group_query, class_query){
  
  all = as.data.frame(count_table) %>%
    filter(group == 'all')
  
  all_class = filter(all, class == class_query)
  
  of_group = as.data.frame(count_table) %>%
    filter(group == group_query)
  
  up_of_group = of_group %>% filter(class == class_query)
  
  N = sum(all$sum)
  k = all_class$sum
  M = sum(of_group$sum)
  n = up_of_group$sum
  
  enrichment = (n * N) / (k * M)

  p = dhyper(x = c(0:n), m = k, n = (N-k), k = M)
  
  return(enrichment)
  
}

calculate_enrichment_factor(group_query = 'upregulated', class_query = 'LTR')

#################################################################
# Gene set enrichment with 'fgsea'
#################################################################

library(fgsea)

## Building list of TEs ranked by p-value

input = results_df

input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))

input = mutate(input, ranking_value = sign(log2FoldChange) * -log10(padj))

ranks = input$ranking_value
names(ranks) = input$ID
ranks = sort(ranks)

## Building pathways

length = length(unique(input$family))

pathways = vector(mode = 'list', length = length)

for (i in unique(input$family)){

  pathways[[i]] = filter(input, family == i)$ID
  
}

pathways = compact(pathways)

## GSEA

fgseaRes = fgsea(pathways = pathways, stats = ranks)

filter(fgseaRes, padj < 0.05)

enrichment_plot = plotEnrichment(pathways[["ERV1"]], ranks) +
  ggtitle('ERV1') +
  xlab('Rank (by p-value)') +
  ylab('Enrichment score') +

enrichment_plot + 
  theme_bw() + 
  theme(plot.title = element_text(face = 'bold', size = 20),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        axis.title = element_text(size = 14)) 

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_ERV1enrichment.png", 
       width = 20, height = 15, units = "cm")