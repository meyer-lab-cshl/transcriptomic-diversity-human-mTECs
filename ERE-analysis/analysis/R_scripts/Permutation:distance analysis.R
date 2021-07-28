library(tidyverse)
library(GenomicRanges)
library(regioneR)
library(pheatmap)

GRanges_gene_extended = loadRDS(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_gene_extended')
GRanges_TE_start = loadRDS(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_TE_start')

#################################################################
# Permutation/distance analysis
#################################################################

## Overlap with GOIs

group_A = list('AIRE_genes' = subset(GRanges_gene_extended, Geneid %in% AIRE_genes),
                 'FEZF2_genes' = subset(GRanges_gene_extended, Geneid %in% FEZF2_genes),
                 'Housekeeping_genes' = subset(GRanges_gene_extended, Geneid %in% housekeeping_genes), 
                 'other_genes' = subset(GRanges_gene_extended, !(Geneid %in% AIRE_genes) & !(Geneid %in% FEZF2_genes) & !(Geneid %in% housekeeping_genes)))

group_B = list('up_EREs' = subset(GRanges_ERE_start, locus %in% up_EREs),
                         'unchanged_EREs' = subset(GRanges_ERE_start, !(locus %in% up_EREs) & !(locus %in% down_EREs)),
                         'down_EREs' = subset(GRanges_ERE_start, locus %in% down_EREs))

perm_test_output_A = run_perm_test(group_A = group_A,
                                   group_B = group_B,
                                   universe = GRanges_gene_extended,
                                   mode = 'overlap',
                                   iterations = 100)

## Overlap with DEGs

group_A = list('up_genes' = subset(GRanges_gene_extended, Geneid %in% up_genes),
               'unchanged_genes' = subset(GRanges_gene_extended, !(Geneid %in% up_genes) & !(Geneid %in% down_genes)),
               'down_genes' = subset(GRanges_gene_extended, Geneid %in% down_genes))

group_B = list('up_EREs' = subset(GRanges_ERE_start, locus %in% up_EREs),
               'unchanged_EREs' = subset(GRanges_ERE_start, !(locus %in% up_EREs) & !(locus %in% down_EREs)),
               'down_EREs' = subset(GRanges_ERE_start, locus %in% down_EREs))

perm_test_output_B = run_perm_test(group_A = group_A,
                                   group_B = group_B,
                                   universe = GRanges_gene_extended,
                                   mode = 'overlap',
                                   iterations = 100)

## Distance from DEGs

group_A = list('up_genes' = subset(GRanges_gene_extended, Geneid %in% up_genes),
               'down_genes' = subset(GRanges_gene_extended, Geneid %in% down_genes))

group_B = list('up_EREs' = subset(GRanges_ERE_start, locus %in% up_EREs),
               'down_EREs' = subset(GRanges_ERE_start, locus %in% down_EREs))

perm_test_output_B = run_perm_test(group_A = group_A,
                                   group_B = group_B,
                                   universe = GRanges_gene_extended,
                                   mode = 'distance',
                                   iterations = 1000)

## Distance from GOIs

group_A = list('AIRE_genes' = subset(GRanges_gene_extended, Geneid %in% AIRE_genes),
               'FEZF2_genes' = subset(GRanges_gene_extended, Geneid %in% FEZF2_genes),
               'Housekeeping_genes' = subset(GRanges_gene_extended, Geneid %in% housekeeping_genes))

group_B = list('up_EREs' = subset(GRanges_ERE_start, locus %in% up_EREs),
               'down_EREs' = subset(GRanges_ERE_start, locus %in% down_EREs))

perm_test_output_A = run_perm_test(group_A = group_A,
                                   group_B = group_B,
                                   universe = GRanges_gene_extended,
                                   mode = 'distance',
                                   iterations = 1000)


## Heatmap

input = perm_test_output_A

labels = input$p_value

for (i in 1:length(labels)){
  
  old_value = labels[i]
  p = input[[2]][i]
  
  significance_threshold = 0.05
  
  if (p < significance_threshold){
    
    p_value = as.character(signif(p, digits = 2))
    annotation = glue('p = {p_value}')
    
  }
  
  else{
    
    annotation = 'ns'
    
  }
  
  new_value = annotation
  labels[i] = new_value
  
}

my_heatmap = pheatmap(mat = -input[[1]],
                      color = colorRampPalette(rev(c("#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB")))(100),
                      cluster_rows=FALSE,
                      show_rownames=TRUE, 
                      cluster_cols=FALSE,
                      fontsize_number = 15,
                      fontsize = 15,
                      display_numbers = labels,
                      number_color = 'black',
                      border_color = 'black',
                      angle_col = 0)

save_pheatmap_png(my_heatmap, "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/GOIs_EREs_distance.png")