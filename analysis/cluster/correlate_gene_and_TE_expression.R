.libPaths('/grid/meyer/home/mpeacey/R/x86_64-pc-linux-gnu-library/4.0/')

GRanges_gene = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/GRanges_TE.rds")
GRanges_TE = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/GRanges_TE.rds")

correlate_fold_change = function(query, subject){
  
  gene_fold_change = rep(NA, length(query))
  TE_fold_change = rep(NA, length(query))
  
  for (gene in c(1:length(query))){
    print(gene)
    gene_fold_change[gene] = query[gene]$log2FoldChange[1]
    
    overlapping_TEs = findOverlaps(query = query[gene],
                                   subject = subject)
    
    overlapping_TEs = as.data.frame(overlapping_TEs)
    hit_indices = overlapping_TEs$subjectHits
    
    fold_changes = rep(NA, length(hit_indices))
    index_number = 1
    
    for (index in hit_indices){
      
      fold_changes[index_number] = 2 ** subject[index]$log2FoldChange[1]
      index_number = index_number + 1
      
    }
    
    TE_fold_change[gene] = log2(mean(fold_changes))
    
  }
  
  output = data.frame(gene_fold_change = gene_fold_change, TE_fold_change = TE_fold_change)
  
  return(output)
  
}

output = correlate_fold_change(query = GRanges_gene, subject = GRanges_TE)

saveRDS(output, "/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/correlate_gene_and_TE_expression_output.rds")


