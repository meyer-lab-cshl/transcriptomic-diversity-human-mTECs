library(tidyverse)
library(GenomicRanges)

input = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_ERE_start')
GRanges_gene_extended = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_gene_extended')
up_genes = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/up_genes')
down_genes = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/down_genes')
diff_genes = append(up_genes, down_genes)

overlaps = as.data.frame(findOverlaps(query = input,
                                      subject = GRanges_gene_extended))

report = vector(length = length(unique(overlaps$queryHits)))

for (entry in 1:length(input)){
  
  if(!(entry %in% overlaps$queryHits)){
    
    report[entry] = 'none'
    
  }
  
  else{
    
    gene_hits = subset(overlaps, queryHits == unique(overlaps$queryHits)[entry])$subjectHits
    
    up = GRanges_gene_extended[gene_hits, ]$Geneid %in% up_genes
    down = GRanges_gene_extended[gene_hits, ]$Geneid %in% down_genes
    unchanged = !(GRanges_gene_extended[gene_hits, ]$Geneid %in% diff_genes)
    
    expression = list('unchanged' = length(unchanged[unchanged == T]),
                      'up' = length(up[up == T]), 
                      'down' = length(down[down == T]))
    
    if(which.max(expression) == 1){
      
      report[entry] = 'unchanged'
      
    }
    
    if(which.max(expression) == 2){
      
      report[entry] = 'up'
      
    }
    
    if(which.max(expression) == 3){
      
      report[entry] = 'down'
      
    }
    
    if(which.max(expression) != 1 &which.max(expression) != 2 & which.max(expression) != 3){
      
      report[entry] = 'other'
      
    }
    
  }
  
  print(entry/length(input) * 100)
  
}

input$overlap_expression = report
saveRDS(input, file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/overlap_annotated_GRanges_gene_extended')