.libPaths = '/grid/meyer/home/mpeacey/R/x86_64-pc-linux-gnu-library/4.0/'

library(GenomicRanges)
library(regioneR)

## Import variables

gene_groups = readRDS("~/TE_thymus/analysis/cluster/gene_groups.rds")
TE_groups = readRDS("~/TE_thymus/analysis/cluster/TE_groups.rds")

##########################

run_perm_test = function(gene_groups, TE_groups, mode = 'distance'){
  
  output = matrix(, nrow = length(gene_groups), ncol = length(TE_groups))
  rownames(output) = c('gene_up', 'gene_down')
  colnames(output) = c('TE_up', 'TE_down')
  
  row_number = 1
  number_of_tests = ncol(output) + nrow(output)
  
  for (gene_group in gene_groups){
    
    column_number = 1
    
    for (TE_group in TE_groups){
      
      print(glue('Starting column {column_number}, row {row_number}'))
      
      if (mode == 'distance'){
        
        pt = permTest(A = TE_group, 
                      B = gene_group, 
                      ntimes = 100,
                      randomize.function = resampleRegions,
                      universe = GRanges_TE,
                      evaluate.function = meanDistance,
                      alternative = 'less')
        
        p_value = pt$meanDistance[[1]]
        Z_score = pt$meanDistance[[6]]
        
      }
      
      if (mode == 'overlap'){
        
        pt = permTest(A = TE_group, 
                      B = gene_group, 
                      ntimes = 100,
                      randomize.function = resampleRegions,
                      universe = GRanges_TE,
                      evaluate.function = numOverlaps,
                      alternative = 'greater')
        
        p_value = pt$numOverlaps[[1]]
        Z_score = -pt$numOverlaps[[6]]
        
      }
      
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
  
  return(output)
  
}

## Run

output = run_perm_test(gene_groups, TE_groups, mode = 'overlap')

# Export output

saveRDS(output, "~/TE_thymus/analysis/cluster/position_enrichment_permutation_output.rds")