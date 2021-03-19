.libPaths('/grid/meyer/home/mpeacey/R/x86_64-pc-linux-gnu-library/4.0/')

library(GenomicRanges)

GRanges_gene = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/GRanges_gene.rds")
GRanges_TE = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/GRanges_TE.rds")

correlate_fold_change = function(query, subject){
  
  query_fold_change = rep(NA, length(query))
  subject_fold_change = rep(NA, length(query))
  
  for (i in 1:length(query)){
    print(i)
    query_fold_change[i] = query[i]$log2FoldChange[1]
    
    overlapping_subjects = findOverlaps(query = query[i],
                                        subject = subject)
    
    overlapping_subjects = as.data.frame(overlapping_subjects)
    hit_indices = overlapping_subjects$subjectHits
    
    fold_changes = rep(NA, length(hit_indices))
    index_number = 1
    
    for (index in hit_indices){
      
      fold_changes[index_number] = 2 ** subject[index]$log2FoldChange[1]
      index_number = index_number + 1
      
    }
    
    subject_fold_change[i] = log2(mean(fold_changes))
    
  }
  
  output = data.frame(query_fold_change = query_fold_change, subject_fold_change = subject_fold_change)
  
  return(output)
  
}

output = correlate_fold_change(query = GRanges_gene, subject = GRanges_TE)

saveRDS(output, "/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/correlate_gene_and_TE_expression_output.rds")


