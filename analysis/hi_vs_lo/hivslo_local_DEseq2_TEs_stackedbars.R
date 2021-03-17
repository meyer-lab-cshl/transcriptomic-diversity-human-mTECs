library(ggplot2)
library(DESeq2)
library(dplyr)
library(tidyr)

#################################################################
# Stacked bar chart: class/family frequency
#################################################################

build_count_table = function(tool, dds, results_df, group, mode){
  
  normalized_counts = as.data.frame(counts(dds, normalized = TRUE))
  
  normalized_counts = cbind(ID = rownames(normalized_counts), normalized_counts)
  
  ## ID column is split into constituent parts, which depend on the tool used.
  
  if (tool == 'TE_transcripts'){
    
    normalized_counts = separate(data = normalized_counts, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
    
  }
  
  if (tool == 'TE_local'){
    
    normalized_counts = separate(data = normalized_counts, col = 'ID', into = c('element', 'gene', 'family', 'class'), sep = ':')
    
  }
  
  ## Row names are reset to ID
  
  normalized_counts = cbind(ID = rownames(normalized_counts), normalized_counts)
  
  ## Input is determined by subsetting normalized_counts by logical conditions specified in 'group'
  
  for (i in group){
    
    if (i == 'all'){
      
      input = normalized_counts
      
    }
    
    if (i == 'diff_regulated'){
      
      sigGenes = results_df[results_df$significant == TRUE,]$ID
      input = normalized_counts[rownames(normalized_counts) %in% sigGenes,]
      
    }
    
    if (i == 'upregulated'){
      
      upGenes = results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange > 0), ]$ID
      input = normalized_counts[rownames(normalized_counts) %in% upGenes,]
      
    }
    
    if (i == 'downregulated'){
      
      downGenes = results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange < 0), ]$ID
      input = normalized_counts[rownames(normalized_counts) %in% downGenes,]
      
    }
    
    if (i == 'diff_regulated_both'){
      
      subset = results_df[results_df$significant == TRUE,]$ID
      input = normalized_counts[rownames(normalized_counts) %in% subset,]
    
      subset_2 = results_df_transcripts[results_df_transcripts$significant == TRUE,]$gene
      input = filter(input, gene %in% subset_2)

    }
    
    input_HI = select(input, -c('214_LO', '221_LO', '226_LO'))
    input_HI$mean = rowMeans(input_HI[6:8])
    input_HI = select(input_HI, -c('214_HI', '221_HI', '226_HI'))
  
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
  
  output = output %>%
    group_by(group) %>%
    mutate(percent_counts = sum / sum(sum) * 100)
  
  return(output)
  
}

## Plot stacked bars

tool = 'TE_local'
dds = dds_local_TE
results_df = results_df_local_TE
group = c('all', 'diff_regulated', 'downregulated', 'upregulated')
mode = 'class'

count_table = build_count_table(tool, dds, results_df, group, mode)

bar_chart = ggplot(count_table, aes(x = group, y = percent_counts, fill = class)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  scale_fill_brewer(palette = "Set1") +
  xlab('') +
  ylab('Fraction of normalized reads') +
  ggtitle('All TE classes', 'TE local')

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

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/TE_local/hi_vs_lo_TEs_classbreakdown.png", 
       width = 20, height = 13, units = "cm")
