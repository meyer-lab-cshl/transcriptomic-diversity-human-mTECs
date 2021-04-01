library(pheatmap)
library(dplyr)

#################################################################
# Heatmap of read counts
#################################################################

## Functions

generate_heatmap_matrix = function(element_mode, tissue_collapse, filter_mode, number_of_elements = 50){
  
  if (element_mode == 'gene'){
    
    input = vs_dds_transcripts_gene
    
  }
  
  if (element_mode == 'TE'){
    
    input = vs_dds_transcripts_TE
    
  }
  
  if (tissue_collapse == T){
    
    input = collapse_tissue_replicates(input)
    
  }
  
  if (tissue_collapse == F){
    
    input = as.data.frame(assay(input))
    
  }
  
  if (filter_mode == 'none'){
    
    matrix = input
    
  }
  
  if (filter_mode == 'variance'){
    
    topVarianceGenes = head(order(rowVars(input), decreasing=T), number_of_elements)
    matrix = input[topVarianceGenes, ]
    
  }
  
  if (filter_mode == 'sig_diff'){
    
    if (element_mode == 'gene'){

      matrix = filter(input, rownames(input) %in% results_df_transcripts_gene_sigdiff$ID)
      
    }
  
    if (element_mode == 'TE'){
      
      matrix = filter(input, rownames(input) %in% results_df_transcripts_TE_sigdiff$ID)
      
    }  
  
  }
  
  return(matrix)
  
}

## LTRs only

matrix = vs_dds_transcripts_TE_collapsed

matrix$ID = rownames(matrix)
matrix = separate(data = matrix, 
                  col = 'ID', 
                  into = c('gene', 'family', 'class'), 
                  sep = ':',
                  remove = T)
matrix = filter(matrix, class == 'LTR')
matrix = select(matrix, -c('gene', 'family', 'class'))

## Plot

matrix = generate_heatmap_matrix(element_mode = 'TE',
                                 tissue_collapse = F,
                                 filter_mode = 'sig_diff')

row_annotation = data.frame(ID = rownames(matrix))
rownames(row_annotation) = row_annotation$ID
row_annotation = separate(data = row_annotation, 
                          col = 'ID', 
                          into = c('gene', 'family', 'class'), 
                          sep = ':',
                          remove = T)
row_annotation = mutate(row_annotation, class = sub("\\?", "", class))
row_annotation = select(row_annotation, class)

my_heatmap = pheatmap(matrix, 
                      cluster_rows=T,
                      show_rownames=F, 
                      cluster_cols=T,
                      scale = 'row',
                      annotation_row = row_annotation)

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "/Users/mpeacey/TE_thymus/analysis/Plots/21-04-01/heatmap_sig_diff_hi_and_lo.png")

