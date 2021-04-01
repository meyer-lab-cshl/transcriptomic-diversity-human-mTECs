library(pheatmap)

#################################################################
# Heatmap 
#################################################################

## Rows ordered by fold change w/out clustering

select <- order(sig_results_df$abs_log2FoldChange, decreasing=TRUE)
pheatmap(sig_transformed_counts[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)

## 1

topVarianceGenes <- head(order(rowVars(assay(vs_dds_transcripts_TE)), decreasing=T),50)

matrix = assay(vs_dds_transcripts_TE)[topVarianceGenes, ]

## 2 

topVarianceGenes = head(order(apply(vs_dds_transcripts_TE_collapsed,1,var), decreasing=T),100)
matrix = vs_dds_transcripts_TE_collapsed[topVarianceGenes, ]

matrix = vs_dds_transcripts_TE_collapsed

## 3

matrix = filter(vs_dds_transcripts_TE_collapsed, 
                rownames(vs_dds_transcripts_TE_collapsed) %in% results_df_transcripts_TE_sigdiff$ID)

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

save_pheatmap_png(my_heatmap, "/Users/mpeacey/TE_thymus/analysis/Plots/21-04-01/heatmap_all_elements.png")

