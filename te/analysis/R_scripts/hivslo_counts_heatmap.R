library(pheatmap)
library(dplyr)

#################################################################
# Heatmap of read counts
#################################################################

## Functions

collapse_tissue_replicates = function(vs_dds, mode = 'mean'){
  
  if (mode == 'mean'){
    
    summary_func = rowMeans
    
  }
  
  if (mode == 'median'){
    
    summary_func = rowMedians
    
  }
  
  counter = 0
  for (i in unique(vs_dds$tissue)){
    
    entry = vs_dds[ , vs_dds$tissue %in% c(i)]
    
    if (counter == 0){
      
      output = data.frame('placeholder' = summary_func(assay(entry)))
      names(output) = i
      
    }
    
    else{
      
      output[i] = summary_func(assay(entry))
      
    }
    
    counter = counter + 1
    
  }
  
  return(output)
  
}

generate_heatmap_matrix = function(input, element_mode, tissue_collapse, filter_mode, number_of_elements = 50){
  
  if (tissue_collapse == T){
    
    input = collapse_tissue_replicates(input, mode = 'mean')
    
  }
  
  if (tissue_collapse == F){
    
    input = as.data.frame(assay(input))
    
  }
  
  if (filter_mode == 'none'){
    
    matrix = input
    
  }
  
  if (filter_mode == 'variance'){
    
    topVarianceGenes = head(order(apply(input,1,var), decreasing=T),number_of_elements)
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


## Row annotation 

row_annotation = data.frame(ID = rownames(matrix))
rownames(row_annotation) = row_annotation$ID
row_annotation = separate(data = row_annotation, 
                          col = 'ID', 
                          into = c('gene', 'family', 'class'), 
                          sep = ':',
                          remove = T)


row_annotation = separate(data = row_annotation, 
                          col = 'ID', 
                          into = c('element', 'gene', 'family', 'class'), 
                          sep = ':',
                          remove = T)

row_annotation = mutate(row_annotation, class = sub("\\?", "", class))
row_annotation = select(row_annotation, class)

## Column annotation

col_annotation = data.frame(ID = colnames(matrix))
rownames(col_annotation) = col_annotation$ID

col_annotation = separate(data = col_annotation, 
                          col = 'ID', 
                          into = c('patient', 'tissue', 'batch'), 
                          sep = '_',
                          remove = T)

col_annotation = select(col_annotation, tissue)

#################################################################
# Heatmap of read counts
#################################################################

matrix = generate_heatmap_matrix(input = vs_dds_transcripts_TE,
                                 element_mode = 'TE',
                                 tissue_collapse = F,
                                 filter_mode = 'none')

my_heatmap = pheatmap(matrix, 
                      cluster_rows=T,
                      show_rownames=F,
                      show_colnames = T,
                      cluster_cols=T,
                      scale = 'row')

## TEtranscripts: top 200 elements by variance

matrix = generate_heatmap_matrix(input = vs_dds_transcripts_TE,
                                 element_mode = 'TE',
                                 tissue_collapse = F,
                                 filter_mode = 'variance',
                                 number_of_elements = 1000)

gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color(6)

ann_colors = list(tissue=c('ESC' = '#F8766D', 
                           'mTEC-hi' = '#B79F00', 
                           'mTEC-lo' = '#00BA38', 
                           'Muscle' = '#00BFC4', 
                           'Ovary' = '#619CFF', 
                           'Testis' = '#F564E3'))

my_heatmap = pheatmap(matrix, 
                      cluster_rows=T,
                      show_rownames=F,
                      show_colnames = T,
                      cluster_cols=T,
                      scale = 'row',
                      annotation_row = row_annotation),
                      annotation_colors = ann_colors,
                      main = 'TE transcripts: top 200 elements by variance')

## TElocal: top 1000 elements by variance

matrix = generate_heatmap_matrix(input = vs_dds_local_TE,
                                 element_mode = 'TE',
                                 tissue_collapse = F,
                                 filter_mode = 'variance',
                                 number_of_elements = 1000)

gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color(6)

ann_colors = list(tissue=c('ESC' = '#F8766D', 
                           'mTEC-hi' = '#B79F00', 
                           'mTEC-lo' = '#00BA38', 
                           'Muscle' = '#00BFC4', 
                           'Ovary' = '#619CFF', 
                           'Testis' = '#F564E3'))

my_heatmap = pheatmap(matrix, 
                      cluster_rows=T,
                      show_rownames=F,
                      show_colnames = F,
                      cluster_cols=T,
                      scale = 'row',
                      annotation_col = col_annotation,
                      annotation_colors = ann_colors,
                      main = 'TE local: top 1000 elements by variance')

save_pheatmap_png <- function(x, filename, width=1200, height=600, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/heatmap_local.png")

