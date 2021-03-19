.libPaths('/grid/meyer/home/mpeacey/R/x86_64-pc-linux-gnu-library/4.0/')

library(GenomicRanges)

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
  
  output = data.frame('query_fold_change' = query_fold_change, 'subject_fold_change' = subject_fold_change)
  
  return(output)
  
}

GRanges_gene_sigdiff = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/correlation_analysis/GRanges_gene_sigdiff.rds")
GRanges_TE_sigdiff = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/correlation_analysis/GRanges_TE_sigdiff.rds")

output = correlate_fold_change(query = GRanges_gene_sigdiff, subject = GRanges_TE_sigdiff)

saveRDS(output, "/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/correlation_analysis/gene_sigdiff_vs_TE_sigdiff.rds")