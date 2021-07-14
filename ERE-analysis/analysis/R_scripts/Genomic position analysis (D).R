library(DESeq2)
library(tidyverse)
library(GenomicRanges)
library(regioneR)
library(pheatmap)
library(RColorBrewer)
library(glue)

functions_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_functions/"
functions = c('extract_subset', 'differential_expression', 'process_DESeq2_results', 'make_GRanges')

for (i in functions){
  
  load(glue('{functions_directory}{i}'))
  
}

#################################################################
# DESeq2
#################################################################

count_table_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/"
data = read.table(glue('{count_table_directory}TE_local_hi_vs_lo.cntTable'),header=T,row.names=1)

dds_local = differential_expression(data, design=~patient+tissue)

results_local = results(dds_local, 
                        independentFiltering = F,
                        contrast = c('tissue', 'mTEC-hi', 'mTEC-lo'))

results_local_gene = extract_subset(mode = 'gene', input = results_local)
results_local_TE = extract_subset(mode = 'TE', input = results_local)

results_df_local_gene = process_DESeq2_results(results = results_local_gene, mode = 'Gene')
results_df_local_TE = process_DESeq2_results(results = results_local_TE, mode = 'TE_local')

results_df_local_gene_up = filter(results_df_local_gene, (significant == T) & (log2FoldChange > 0))
results_df_local_gene_unchanged = filter(results_df_local_gene, significant == F)
results_df_local_gene_down = filter(results_df_local_gene, (significant == T) & (log2FoldChange < 0))
results_df_local_gene_sigdiff = filter(results_df_local_gene, significant == T)

results_df_local_TE_up = filter(results_df_local_TE, (significant == T) & (log2FoldChange > 0))
results_df_local_TE_unchanged = filter(results_df_local_TE, significant == F)
results_df_local_TE_down = filter(results_df_local_TE, (significant == T) & (log2FoldChange < 0))
results_df_local_TE_sigdiff = filter(results_df_local_TE, significant == T)

#################################################################
# Genomic position analysis (D)
#################################################################

## Generate GRanges objects

GRanges_TE = make_GRanges(mode = 'TE',
                          results_df = results_df_local_TE)
GRanges_TE_start = GRanges_TE
end(GRanges_TE_start) = GenomicRanges::start(GRanges_TE_start) + 100

GRanges_TE_up = make_GRanges(mode = 'TE',
                             results_df = results_df_local_TE_up)
GRanges_TE_up_start = GRanges_TE_up
end(GRanges_TE_up_start) = GenomicRanges::start(GRanges_TE_up) + 100
start(GRanges_TE_up_start) = GenomicRanges::start(GRanges_TE_up) - 100

GRanges_TE_down = make_GRanges(mode = 'TE',
                               results_df = results_df_local_TE_down)
GRanges_TE_down_start = GRanges_TE_down
end(GRanges_TE_down_start) = GenomicRanges::start(GRanges_TE_down_start) + 100

GRanges_TE_unchanged = make_GRanges(mode = 'TE',
                                    results_df = results_df_local_TE_unchanged)
GRanges_TE_unchanged_start = GRanges_TE_unchanged
end(GRanges_TE_unchanged_start) = GenomicRanges::start(GRanges_TE_unchanged_start) + 100

## A: Overlap between genes and TEs

perm_test_output_A = run_perm_test(group_A = list(TE_up = GRanges_TE_up, TE_unchanged = GRanges_TE_unchanged, TE_down = GRanges_TE_down),
                                   group_B = list(gene_up = GRanges_gene_up, gene_unchanged = GRanges_gene_unchanged, gene_down = GRanges_gene_down),
                                   universe = GRanges_TE,
                                   mode = 'overlap',
                                   iterations = 1000)

labels = signif(perm_test_output_A$evaluation, digits = 3)

for (i in 1:length(labels)){
  
  old_value = labels[i]
  p = perm_test_output_A[[2]][i]
  print(p)
  
  significance_threshold = 0.05 / length(labels)
  
  if (p < significance_threshold){
    
    p_value = as.character(signif(p, digits = 1))
    annotation = glue('p = {p_value}')
    
  }
  
  else{
    
    annotation = 'ns'
    
  }
  
  new_value = glue('{old_value}%, {annotation}')
  labels[i] = new_value
  
}

print(labels)

my_heatmap = pheatmap(mat = perm_test_output_A[[1]],
                      color = colorRampPalette(rev(c("#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB")))(100),
                      cluster_rows=FALSE,
                      show_rownames=TRUE, 
                      cluster_cols=FALSE,
                      display_numbers = labels,
                      fontsize_number = 15,
                      fontsize = 15,
                      number_color = 'black',
                      border_color = 'black',
                      angle_col = '0',
                      labels_row = c('  Up', '   - ', '  Down'),
                      labels_col = c('Up', '-', 'Down'))

save_pheatmap_png(my_heatmap, "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/gene_TE_overlap.png")

## Overlap between TSSs and TEs

df_all = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_all.csv') 
GRanges_TSS = makeGRangesFromDataFrame(df_all, keep.extra.columns = T)

df_high = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_high.csv') 
GRanges_TSS_high = makeGRangesFromDataFrame(df_high, keep.extra.columns = T)

df_low = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_low.csv') 
GRanges_TSS_low = makeGRangesFromDataFrame(df_low, keep.extra.columns = T)

output = run_perm_test(group_A = list(TE_up = GRanges_TE_up, TE_down = GRanges_TE_down),
                       group_B = list(TSS_high = GRanges_TSS_high, TSS_low = GRanges_TSS_low),
                       universe = GRanges_TE,
                       mode = 'overlap')
