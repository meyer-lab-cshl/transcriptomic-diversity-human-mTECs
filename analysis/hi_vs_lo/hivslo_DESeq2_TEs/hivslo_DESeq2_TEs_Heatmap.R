#################################################################
# Heatmap version 1 (ggplot)
#################################################################

sig_normalized_counts_long <- melt(sig_normalized_counts, id.vars=c("ID"))

ggplot(sig_normalized_counts_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_raster() + theme(axis.text.x=element_text(angle=65, hjust=1)) 

#################################################################
# Heatmap version 2 (pheatmap)
#################################################################

## Rows ordered by fold change w/out clustering

select <- order(sig_results_df$abs_log2FoldChange, decreasing=TRUE)
pheatmap(sig_transformed_counts[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)

## Rows clustered

my_heatmap = pheatmap(sig_transformed_counts, 
                      cluster_rows=TRUE,
                      show_rownames=TRUE, 
                      cluster_cols=TRUE,
                      cutree_cols = 2)

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_heatmap.png")