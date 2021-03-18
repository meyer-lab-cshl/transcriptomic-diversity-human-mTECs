## Import variables

gene_groups = readRDS("~/TE_thymus/analysis/cluster_R_stuff/gene_groups.rds")
TE_groups = readRDS("~/TE_thymus/analysis/cluster_R_stuff/TE_groups.rds")

##########################

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


# Export output

saveRDS(output, "~/TE_thymus/analysis/cluster_R_stuff/gene_groups.rds")