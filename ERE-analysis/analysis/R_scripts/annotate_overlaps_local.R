library(tidyverse)
library(GenomicRanges)
<<<<<<< HEAD
library(glue)

working_directory = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis'

results_df_local_ERE = readRDS(file = glue('{working_directory}/R_variables/results_df_local_ERE'))
input = readRDS(file = glue('{working_directory}/R_variables/GRanges_ERE_start')) %>%
  subset(locus %in% subset(results_df_local_ERE, significant == T)$locus)
GRanges_gene_extended = readRDS(file = glue('{working_directory}/R_variables/GRanges_gene_extended'))
up_genes = readRDS(file = glue('{working_directory}/R_variables/up_genes'))
down_genes = readRDS(file = glue('{working_directory}/R_variables/down_genes'))
=======

results_df_local_ERE = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/results_df_local_ERE')
input = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_ERE_start') %>%
  subset(locus %in% results_df_local_ERE$locus)
GRanges_gene_extended = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_gene_extended')
up_genes = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/up_genes')
down_genes = readRDS(file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/down_genes')
>>>>>>> 73018c40b61fd01ae27c86a356ad3cdcf97896c3
diff_genes = append(up_genes, down_genes)

overlaps = as.data.frame(findOverlaps(query = input,
                                      subject = GRanges_gene_extended))

report = vector(length = length(unique(overlaps$queryHits)))

for (entry in 1:length(input)){
  
  if(!(entry %in% overlaps$queryHits)){
    
    report[entry] = 'none'
    
  }
  
  else{
    
    gene_hits = subset(overlaps, queryHits == entry)$subjectHits
    
    up = GRanges_gene_extended[gene_hits, ]$Geneid %in% up_genes
    down = GRanges_gene_extended[gene_hits, ]$Geneid %in% down_genes
    unchanged = !(GRanges_gene_extended[gene_hits, ]$Geneid %in% diff_genes)
    
    expression = list('up' = length(up[up == T]), 
                      'down' = length(down[down == T]),
                      'unchanged' = length(unchanged[unchanged == T]))
    
    if(which.max(expression) == 1){
      
      report[entry] = 'up'
      
    }
    
    if(which.max(expression) == 2){
      
      report[entry] = 'down'
      
    }
    
    if(which.max(expression) == 3){
      
      report[entry] = 'unchanged'
      
    }
    
    if(which.max(expression) != 1 &which.max(expression) != 2 & which.max(expression) != 3){
      
      report[entry] = 'other'
      
    }
    
  }
  
  print(entry/length(input) * 100)
  
}

input$overlap_expression = report
<<<<<<< HEAD
saveRDS(input, file = glue('{working_directory}/R_variables/overlap_annotated_GRanges_gene_extended'))
=======
saveRDS(input, file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/overlap_annotated_GRanges_gene_extended')
>>>>>>> 73018c40b61fd01ae27c86a356ad3cdcf97896c3
