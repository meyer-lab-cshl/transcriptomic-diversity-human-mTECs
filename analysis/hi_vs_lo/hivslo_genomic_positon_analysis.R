library(GenomicRanges)
library(regioneR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggplotify)
library(pheatmap)
library(RColorBrewer)

#################################################################
# 
#################################################################

make_GRanges = function(mode, results_df){
  
  if (mode == 'TE'){
    
    annotation = read.table(file = "/Users/mpeacey/TE_thymus/analysis/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt", header = 1)
    annotation = separate(annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':')
    annotation = separate(annotation, start.stop, into = c('start', 'end'), sep = '-')
    annotation = rename(annotation, locus = TE)
    
    df = merge(results_df, annotation, by = 'locus')
    
  }
  
  if (mode == 'gene'){
    
    annotation = read.table(file = '/Users/mpeacey/TE_thymus/analysis/annotation_tables/gencode.v38_gene_annotation_table.txt', header = 1)
    annotation = select(annotation, c('Geneid', 'Chromosome', 'Start', 'End', 'Strand'))
    annotation = rename(annotation, chr = Chromosome, start = Start, end = End, strand = Strand)
    
    df = merge(results_df, annotation, by = 'Geneid')
    
  }
  
  output = makeGRangesFromDataFrame(df, keep.extra.columns = T)
  
  return(output)
  
}

## GRanges


GRanges_gene = make_GRanges(mode = 'gene',
                            results_df = results_df_local_gene)

GRanges_TE = make_GRanges(mode = 'TE',
                          results_df = results_df_local_TE)



GRanges_gene_sigdiff = make_GRanges(mode = 'gene',
                                    results_df = results_df_local_gene_sigdiff)


GRanges_TE_sigdiff = make_GRanges(mode = 'TE',
                                  results_df = results_df_local_TE_sigdiff)

subsetByOverlaps(GRanges_TE_sigdiff, GRanges_gene_sigdiff)

#################################################################
# GRanges
################################################################### 

results_df_local_gene_up = filter(results_df_local_gene, (significant == T) & (log2FoldChange > 0))
results_df_local_gene_unchanged = filter(results_df_local_gene, significant == F)
results_df_local_gene_down = filter(results_df_local_gene, (significant == T) & (log2FoldChange < 0))

results_df_local_TE_up = filter(results_df_local_TE, (significant == T) & (log2FoldChange > 0))
results_df_local_TE_unchanged = filter(results_df_local_TE, significant == F)
results_df_local_TE_down = filter(results_df_local_TE, (significant == T) & (log2FoldChange < 0))

GRanges_gene_up = make_GRanges(mode = 'gene',
                               results_df = results_df_local_gene_up)
GRanges_gene_unchanged = make_GRanges(mode = 'gene',
                               results_df = results_df_local_gene_unchanged)
GRanges_gene_down = make_GRanges(mode = 'gene',
                               results_df = results_df_local_gene_down)

GRanges_TE_up = make_GRanges(mode = 'TE',
                               results_df = results_df_local_TE_up)
GRanges_TE_unchanged = make_GRanges(mode = 'TE',
                                      results_df = results_df_local_TE_unchanged)
GRanges_TE_down = make_GRanges(mode = 'TE',
                                 results_df = results_df_local_TE_down)

#################################################################
# regioneR
################################################################### 

random_RS = resampleRegions(GRanges_TE_sigdiff, universe=GRanges_TE)

## Do my differentially expressed TEs overlap with differentially expressed genes more frequently than expected by chance?

pt = permTest(A = GRanges_TE_sigdiff, 
              B = GRanges_gene_sigdiff, 
              ntimes = 100,
              randomize.function = resampleRegions,
              universe = GRanges_TE,
              evaluate.function = numOverlaps)

## Are differentially expressed TEs closer to differentially expressed genes than expected by chance?

pt = permTest(A = GRanges_TE_sigdiff, 
              B = GRanges_gene_sigdiff, 
              ntimes = 100,
              randomize.function = resampleRegions,
              universe = GRanges_TE,
              evaluate.function = meanDistance)

pt_df = data.frame(permutations = pt$meanDistance[[5]])
              
ggplot() +
  geom_histogram(data = pt_df, aes(x = permutations)) 

## Heatmap

#gene_groups = list(GRanges_gene_up, GRanges_gene_unchanged, GRanges_gene_down)
#TE_groups = list(GRanges_TE_up, GRanges_TE_unchanged, GRanges_TE_down)
gene_groups = list(GRanges_gene_up, GRanges_gene_down)
TE_groups = list(GRanges_TE_up, GRanges_TE_down)

output = matrix(, nrow = length(gene_groups), ncol = length(TE_groups))
rownames(output) = c('gene_up', 'gene_down')
colnames(output) = c('TE_up', 'TE_down')

row_number = 1
number_of_tests = ncol(output) + nrow(output)

for (gene_group in gene_groups){
  
  column_number = 1
  
  for (TE_group in TE_groups){
    
    pt = permTest(A = TE_group, 
                 B = gene_group, 
                 ntimes = 1000,
                 randomize.function = resampleRegions,
                 universe = GRanges_TE,
                 evaluate.function = meanDistance,
                 alternative = 'less')
    
    p_value = pt$meanDistance[[1]]
    Z_score = pt$meanDistance[[6]]
    
    if (p_value < (0.05/number_of_tests)){
      
        output[row_number, column_number] = Z_score
        
      }
  
    else{
      
      output[row_number, column_number] = NA
      
    }
    
    column_number = column_number + 1
    
  }
  
  row_number = row_number + 1
  
}

my_heatmap = pheatmap(mat = output, 
                      cluster_rows=FALSE,
                      show_rownames=TRUE, 
                      cluster_cols=FALSE,
                      color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                "Blues")))(100),
                      display_numbers = T,
                      fontsize_number = 15,
                      number_color = 'black')

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "/Users/mpeacey/TE_thymus/analysis/Plots/TE_local/hi_vs_lo_genomic_poisitons_heatmap.png")