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
  normalized_counts = select(normalized_counts, -c('214_LO', '221_LO', '226_LO'))
  normalized_counts$mean = rowMeans(normalized_counts[2:4])
  normalized_counts = select(normalized_counts, -c('214_HI', '221_HI', '226_HI'))
  
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

count_table = build_count_table(dds = dds_local_TE,
                                results_df = results_df_local_TE, 
                                group = c('all', 'diff_regulated', 'up_regulated', 'down_regulated'), 
                                mode = 'overlap',
                                by = 'locus')

## Plot stacked bars

bar_chart = ggplot(count_table, aes(x = group, y = percent, fill = overlap_status)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  scale_fill_brewer(palette = "Set1") +
  xlab('') +
  ylab('Fraction of normalized reads') +
  ggtitle('TEs in mTEC-HI cells', 'Subset by physical overlap with genes')

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

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/21-03-25/slide_3.png", 
       width = 20, height = 13, units = "cm")
