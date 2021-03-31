library(pheatmap)

#################################################################
# Heatmap version 2 (pheatmap)
#################################################################

## Rows ordered by fold change w/out clustering

select <- order(sig_results_df$abs_log2FoldChange, decreasing=TRUE)
pheatmap(sig_transformed_counts[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)

## Rows clustered

my_heatmap = pheatmap(transformed_counts_TE, 
                      cluster_rows=TRUE,
                      show_rownames=FALSE, 
                      cluster_cols=TRUE)

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_heatmap.png")

#################################################################
# Heatmap version 3
#################################################################

## 1

topVarianceGenes <- head(order(rowVars(assay(vs_dds_transcripts_TE)), decreasing=T),50)

matrix = assay(vs_dds_transcripts_TE)[topVarianceGenes, ]

my_heatmap = pheatmap(matrix, 
                      cluster_rows=TRUE,
                      show_rownames=FALSE, 
                      cluster_cols=TRUE,
                      scale = 'row')

##2 

topVarianceGenes = head(order(apply(vs_dds_transcripts_TE_collapsed,1,var), decreasing=T),100)
matrix = vs_dds_transcripts_TE_collapsed[topVarianceGenes, ]

matrix = vs_dds_transcripts_TE_collapsed

my_heatmap = pheatmap(matrix, 
                      cluster_rows=T,
                      show_rownames=FALSE, 
                      cluster_cols=T,
                      scale = 'row')
