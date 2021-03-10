library(ggplot2)
library(plyr)
library(dplyr)
library(glue)

#################################################################
# Gene set enrichment with 'fgsea'
#################################################################

library(fgsea)

## Building list of TEs ranked by p-value

input = results_df

input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))

input = filter(input, class != 'Unknown')

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

## Plot

enrichment_plot = function(query_pathway){
  
  p_value = filter(fgseaRes, pathway == query_pathway)[, padj]
  p_value = signif(p_value, digits = 3)
  
  NES = filter(fgseaRes, pathway == query_pathway)[, NES]
  NES = signif(NES, digits = 3)
  
  enrichment_plot = plotEnrichment(pathways[[query_pathway]], ranks) +
    ggtitle(query_pathway) +
    xlab('Rank (by p-value)') +
    ylab('Enrichment score') +
    annotate('label', x = Inf, y = Inf, vjust = 3, hjust = 1.3,
             label = glue('Adjusted p-value: {p_value}'))
  
  return(enrichment_plot + 
    theme_bw() + 
    theme(plot.title = element_text(face = 'bold', size = 20),
          panel.border = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(size = 0.8),
          axis.title = element_text(size = 14)))
  
  ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_{query_pathway}enrichment.png", 
         width = 20, height = 15, units = "cm")
}

enrichment_plot('Satellite')

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_ERV1enrichment.png", 
       width = 20, height = 15, units = "cm")

#################################################################
# Enrichment (depricated)
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
