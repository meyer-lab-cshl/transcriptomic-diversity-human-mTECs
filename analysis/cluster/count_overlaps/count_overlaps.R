.libPaths('/grid/meyer/home/mpeacey/R/x86_64-pc-linux-gnu-library/4.0/')

library(GenomicRanges)

GRanges_gene_up = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/objects/GRanges_gene_up.rds")
GRanges_gene_down = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/objects/GRanges_gene_down.rds")
GRanges_TE = readRDS("/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/objects/GRanges_TE.rds")

overlap_with_up_gene = rep(NA, length(GRanges_TE))
overlap_with_down_gene = rep(NA, length(GRanges_TE))

for (i in 1:length(GRanges_TE)){
  
  print(i)
  up_hit = countOverlaps(query = GRanges_TE[i], subject = GRanges_gene_up)
  if (length(up_hit) > 0){
    
    overlap_with_up_gene[i] = T
    
  }
  
  else{
    
    overlap_with_up_gene[i] = F
    
  }
  
  down_hit = countOverlaps(query = GRanges_TE[i], subject = GRanges_gene_down)
  if (length(down_hit) > 0){
    
    overlap_with_down_gene[i] = T
    
  }
  
  else{
    
    overlap_with_down_gene[i] = F
    
  }
  
}

saveRDS(overlap_with_up_gene, "/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/objects/overlap_with_up_gene.rds")
saveRDS(overlap_with_down_gene, "/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/objects/overlap_with_down_gene.rds")
