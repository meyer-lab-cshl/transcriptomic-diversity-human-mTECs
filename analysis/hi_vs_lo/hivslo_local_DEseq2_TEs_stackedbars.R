library(ggplot2)
library(DESeq2)
library(dplyr)
library(tidyr)

#################################################################
# Stacked bar chart: class/family frequency
#################################################################

build_count_table = function(dds, results_df, group, mode, by){
  
  normalized_counts = as.data.frame(counts(dds, normalized = TRUE))
  
  normalized_counts = cbind(ID = rownames(normalized_counts), normalized_counts)
  
  ## Merge with results_df, determine mean and remove unnecessary columns
  
  normalized_counts = merge(normalized_counts, results_df, by = 'ID')
  normalized_counts = select(normalized_counts, -c('pt214_mTEC-lo_new', 'pt221_mTEC-lo_new', 'pt226_mTEC-lo_new'))
  normalized_counts$mean = rowMeans(normalized_counts[2:4])
  normalized_counts = select(normalized_counts, -c('pt214_mTEC-hi_new', 'pt221_mTEC-hi_new', 'pt226_mTEC-hi_new'))
  
  normalized_counts = mutate(normalized_counts, class = sub("\\?", "", class))
  
  ## Input is determined by subsetting normalized_counts by logical conditions specified in 'group'
  
  for (i in group){
    
    if (i == 'all'){
      
      input = normalized_counts
      
    }
    
    if (i == 'diff_regulated'){
      
      sigGenes = results_df[results_df$significant == TRUE,]$ID
      input = normalized_counts[normalized_counts$ID %in% sigGenes,]
      
    }
    
    if (i == 'up_regulated'){
      
      upGenes = results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange > 0), ]$ID
      input =normalized_counts[normalized_counts$ID %in% upGenes,]
      
    }
    
    if (i == 'down_regulated'){
      
      downGenes = results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange < 0), ]$ID
      input = normalized_counts[normalized_counts$ID %in% downGenes,]
      
    }
    
    if (i == 'diff_regulated_both'){
      
      subset = results_df[results_df$significant == TRUE,]$ID
      input = normalized_counts[rownames(normalized_counts) %in% subset,]
    
      subset_2 = results_df_transcripts[results_df_transcripts$significant == TRUE,]$gene
      input = filter(input, gene %in% subset_2)

    }
    
    ########
    
    if (mode == 'class'){
      
      if (by == 'normalized_reads'){
        
        input = group_by(input, class) %>% summarize(summary = sum(mean))
        
      }
      
      if (by == 'locus'){
        
        input = group_by(input, class) %>% summarize(summary = n())
        
      }
      
    }
      
    if (mode == 'LTR_family'){
      
      if (by == 'normalized_reads'){
        
        input = input %>% 
          filter(class == 'LTR') %>%
          group_by(family) %>% 
          summarize(summary = sum(mean))
        
      }
      
      if (by == 'locus'){
        
        input = input %>% 
          filter(class == 'LTR') %>%
          group_by(family) %>% 
          summarize(summary = n())
      
      }
      
    }
    
    if (mode == 'overlap'){
      
      if (by == 'normalized_reads'){
        
        input = group_by(input, overlap_status) %>% summarize(summary = sum(mean))
        
      }
      
      if (by == 'locus'){
        
        input = group_by(input, overlap_status) %>% summarize(summary = n())
        
      }
      
    }
    
    
    
    
    ######
    
    input = as.data.frame(input)
    input$group = i
    
    if (match(i, group) == 1){
      
      output = input
      
    }
    
    if (match(i, group) != 1) {
      
      output = bind_rows(output, input)
      
    }  

  }
  
  output = output %>%
    group_by(group) %>%
    mutate(percent = summary / sum(summary) * 100)
  
  return(output)
  
}


build_count_table = function(group){
  
  ## Input is determined by subsetting normalized_counts by logical conditions specified in 'group'
  
  for (i in group){
    
    if (i == 'all'){
      
      input = as.data.frame(GRanges_TE_annotated)
      
    }
    
    if (i == 'diff_regulated'){
      
      input = as.data.frame(GRanges_TE_annotated) %>%
        filter(significant == TRUE)
      
    }
    
    if (i == 'up_regulated'){
      
      input = as.data.frame(GRanges_TE_annotated) %>%
        filter((significant == TRUE) & (log2FoldChange > 0))
      
    }
    
    input = group_by(input, Class) %>% summarize(summary = n())
    
    ######
    
    input = as.data.frame(input)
    input$group = i
    
    if (match(i, group) == 1){
      
      output = input
      
    }
    
    if (match(i, group) != 1) {
      
      output = bind_rows(output, input)
      
    }  
    
  }
  
  output = output %>%
    group_by(group) %>%
    mutate(percent = summary / sum(summary) * 100)
  
  return(output)
  
}

count_table = build_count_table(dds_local_TE, 
                                results_df_local_TE, 
                                group = c('all', 'down_regulated', 'up_regulated'),
                                mode = 'class',
                                by = 'normalized_reads')

## Plot stacked bars

bar_chart = ggplot(count_table, aes(x = group, y = percent, fill = class)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  scale_fill_brewer(palette = 'Set1') +
  xlab('mTEC-HI subset') +
  ylab('Fraction of normalized reads') +
  labs(fill= "Class") +
  scale_x_discrete(labels = c('All', 'Down', 'Up'))

  
  #ggtitle('TEs in mTEC-HI cells', 'Subset by physical overlap with genes')

bar_chart + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                               plot.subtitle = element_text(size = 14),
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_line(color = 'grey'),
                               panel.grid.minor.y = element_blank(),
                               axis.text.x = element_text(size = 13, margin = margin(t = 6)),
                               axis.text.y = element_text(size = 14),
                               axis.title.y = element_text(size = 14),
                               axis.title.x = element_text(size = 14, margin = margin(t = 6)),
                               axis.line = element_line(size = 0.8),
                               panel.border = element_blank(),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 14))

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/local_class_distribution.png", 
       width = 8, height = 4.5, units = "in")



## Chi-square

contingency = pivot_wider(count_table, id_cols = class, names_from = group, values_from = percent) %>% 
  mutate_if(is.numeric, list(~replace_na(., 0))) %>% 
  as.data.frame() %>%
  pivot_longer(cols = c(all, down_regulated, up_regulated), names_to = 'group', values_to = 'count')

contingency_table = xtabs(count ~ group+class, data=contingency)

vcd::mosaic(~group+class, data = contingency_table, direction = c('v', 'h'), shade = T)

## Chi-square v2

contingency = pivot_wider(count_table, id_cols = class, names_from = group, values_from = percent) %>% 
  mutate_if(is.numeric, list(~replace_na(., 0))) %>% 
  as.data.frame()

rownames(contingency) = contingency$class
contingency = dplyr::select(contingency, -class)

chisq.test(contingency)

chisq.posthoc.test::chisq.posthoc.test(contingency)


