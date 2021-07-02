.libPaths('/grid/meyer/home/mpeacey/R/x86_64-pc-linux-gnu-library/4.0/')

library(DESeq2)
library(dplyr)
library(tidyr)

#################################################################
# Functions
#################################################################

differential_expression = function(raw_count_table, min_reads = 10){
  
  ## Define sampleInfo
  
  ID = colnames(raw_count_table)
  sampleInfo = data.frame(ID,row.names=colnames(raw_count_table))
  sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'group'), sep = '_'))
  sampleInfo$patient = factor(sampleInfo$patient)
  
  ## Construct DESeq dataset object
  
  dds <- DESeqDataSetFromMatrix(countData = raw_count_table, colData = sampleInfo, design = ~patient + group)
  dds$group = relevel(dds$group,ref="LO")
  
  ## Run differential expression analysis
  
  dds = DESeq(dds)
  
  ## Pre-filter to remove rows with less than 'min_reads' total (default 10)
  
  dds = dds[ rowSums(counts(dds)) >= min_reads, ]
  
  return(dds)
  
}

## Split into gene and TE

extract_from_DESeq2 = function(mode, input){
  
  if (mode == 'gene'){
    
    output = input[grep("^ENSG",rownames(input)),]
    
  }
  
  if (mode == 'TE'){
    
    output = input[grepl("^(?!ENSG).*$",rownames(input), perl = TRUE),]
    
  }
  
  return(output)
  
}

## Process results

process_DESeq2_results = function(results,
                                  mode,
                                  p_value_cutoff = 0.1, 
                                  log_fold_change_cutoff = 0.58){
  
  ## Add statistical and biological significance labels
  
  results_df = as.data.frame(results)
  
  results_df = mutate(results_df, significant = case_when(padj < p_value_cutoff ~ TRUE, 
                                                          padj >= p_value_cutoff ~ FALSE))
  
  results_df = mutate(results_df, FC_significant = case_when(abs(log2FoldChange) > log_fold_change_cutoff ~ TRUE, 
                                                             abs(log2FoldChange) <= log_fold_change_cutoff ~ FALSE))
  
  results_df = mutate(results_df, overall_significant = case_when((significant == TRUE) & (FC_significant == TRUE) ~ TRUE, 
                                                                  (significant == FALSE) | (FC_significant == FALSE) ~ FALSE))
  
  ## Add separate ID columns
  
  if (mode == 'TE_transcripts'){
    
    results_df = cbind(ID = rownames(results_df), results_df)
    results_df = separate(data = results_df, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
    results_df = cbind(ID = rownames(results_df), results_df)
    
  }
  
  if (mode == 'TE_local'){
    
    results_df = cbind(ID = rownames(results_df), results_df)
    results_df = separate(data = results_df, col = 'ID', into = c('locus', 'gene', 'family', 'class'), sep = ':')
    results_df = cbind(ID = rownames(results_df), results_df)
    
  }
  
  if (mode == 'Gene'){
    
    results_df = cbind(Geneid = rownames(results_df), results_df)
    
  }
  
  return(results_df)
  
}


#################################################################
# Analysis
#################################################################

data = read.table("/grid/meyer/home/mpeacey/TE_thymus/analysis/count_tables/TE_local_hi_vs_lo.cntTable",
                  header=T,row.names=1)
colnames(data) = c('214_HI', '214_LO', '221_HI', '221_LO', '226_HI', '226_LO')

min_reads = c(2, 5, 10, 20)
independent_filtering = c(T, F)

gene = matrix(, nrow = length(independent_filtering), ncol = length(min_reads))
rownames(gene) = c('IF_on', 'IF_off')
colnames(gene) = c('2', '5', '10', '20')

TE = matrix(, nrow = length(independent_filtering), ncol = length(min_reads))
rownames(TE) = c('IF_on', 'IF_off')
colnames(TE) = c('2', '5', '10', '20')

row_number = 1

for (i in independent_filtering){
  
  column_number = 1
  
  for (i_2 in min_reads){
    
    dds_local = differential_expression(data, min_reads = i_2)
    
    results_local = results(dds_local, independentFiltering = i)
    
    results_local_gene = extract_from_DESeq2(mode = 'gene', input = results_local)
    results_local_TE = extract_from_DESeq2(mode = 'TE', input = results_local)
    
    results_df_local_gene = process_DESeq2_results(results = results_local_gene, mode = 'Gene')
    results_df_local_TE = process_DESeq2_results(results = results_local_TE, mode = 'TE_local')
    
    results_df_local_gene_sigdiff = filter(results_df_local_gene, significant == T)
    
    results_df_local_TE_sigdiff = filter(results_df_local_TE, significant == T)
    
    gene[row_number, column_number] = nrow(results_df_local_gene_sigdiff)
    TE[row_number, column_number] = nrow(results_df_local_TE_sigdiff)
    
    column_number = column_number + 1
    
  }
  
  row_number = row_number + 1
  
}

saveRDS(gene, "/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/objects/sig_diff_gene_count.rds")
saveRDS(TE, "/grid/meyer/home/mpeacey/TE_thymus/analysis/cluster/objects/sig_diff_TE_count.rds")