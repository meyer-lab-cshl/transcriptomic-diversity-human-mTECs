#################################################################
# Stacked bar chart: class/family frequency
#################################################################

build_count_table = function(group, mode){
  
  normalized_counts = as.data.frame(counts(dds, normalized = TRUE))
  
  for (i in group){
    
    if (i == 'all'){
      
      input = normalized_counts
      
    }
    
    if (i == 'diff_regulated'){
      
      input = normalized_counts[rownames(normalized_counts) %in% sigGenes,]
      
    }
    
    if (i == 'upregulated'){
      
      input = normalized_counts[rownames(normalized_counts) %in% upGenes,]
      
    }
    
    if (i == 'downregulated'){
      
      input = normalized_counts[rownames(normalized_counts) %in% downGenes,]
      
    }
    
    input_HI = select(input, c('214_HI', '221_HI', '226_HI'))
    input_HI$mean = rowMeans(input_HI)
    input_HI = select(input_HI, mean)
    
    input_HI = cbind(ID = rownames(input_HI), input_HI)
    input_HI = separate(data = input_HI, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
    input_HI = cbind(ID = rownames(input_HI), input_HI)
    
    input_HI = mutate(input_HI, class = sub("\\?", "", class))
    
    if (mode == 'class'){
      
      input_HI = group_by(input_HI, class) %>% summarize(sum = sum(mean))
      input_HI = as.data.frame(input_HI)
      input_HI$group = i
      
      if (match(i, group) == 1){
        
        output = input_HI
        
      }
      
      else{
        
        output = bind_rows(output, input_HI)
        
      }
      
      if (mode == 'LTR_family'){
        
        input_HI = input_HI %>% 
          filter(class == 'LTR') %>%
          group_by(family) %>% 
          summarize(sum = sum(mean))
        input_HI = as.data.frame(input_HI)
        input_HI$group = i
        
        print(input_HI)
        
        if (match(i, group) == 1){
          
          output = input_HI
          
        }
        
        else{
          
          output = bind_rows(output, input_HI)
          
        }
        
      }  
      
    }
    
  }
  
  output = output %>%
    group_by(group) %>%
    mutate(percent_counts = sum / sum(sum) * 100)
  
  return(output)
  
}

normalized_counts_HI_LTR_family = normalized_counts_HI %>% 
  filter(class == 'LTR') %>%
  group_by(family) %>% 
  summarize(sum = sum(mean))
normalized_counts_HI_LTR_family = as.data.frame(normalized_counts_HI_LTR_family)
normalized_counts_HI_LTR_family$group = 'All'


## Plot stacked bars

group = c('all', 'diff_regulated', 'upregulated', 'downregulated')
mode = 'class'

count_table = build_count_table(group, mode)

bar_chart = ggplot(count_table, aes(x = group, y = percent_counts, fill = class)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  scale_fill_brewer(palette = "Set1") +
  xlab('') +
  ylab('Fraction of normalized reads')

bar_chart + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
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

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_classbreakdown.png", 
       width = 20, height = 15, units = "cm")

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
  
  return(enrichment)
  
}

calculate_enrichment_factor(group_query = 'upregulated', class_query = 'LTR')