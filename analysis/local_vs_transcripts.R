#################################################################
# Count table produced by TE_local
#################################################################

data = read.table("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo_local/count_table/TE_local_hi_vs_lo.cntTable",
                  header=T,row.names=1)
colnames(data) = c('214_HI', '214_LO', '221_HI', '221_LO', '226_HI', '226_LO')

TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

min_read = 1
input = TE_data[apply(TE_data,1,function(x){max(x)}) > min_read,]

input = cbind(ID = rownames(input), input)
input = separate(data = input, col = 'ID', into = c('element', 'gene', 'family', 'class'), sep = ':')
input = cbind(ID = rownames(input), input)

index = select(input, c('ID', 'gene', 'family', 'class'))
index = unique(index[, 2:4])

data_local = aggregate(cbind(input$'214_HI', input$'214_LO', 
                input$'221_HI', input$'221_LO', 
                input$'226_HI', input$'226_LO'), 
                by = list(gene = input$gene), FUN = sum)

data_local = merge(data_local, index, by = 'gene')
data_local = rename(data_local, '214_HI' = V1, '214_LO' = V2,
              '221_HI' = V3, '221_LO' = V4,
              '226_HI' = V5, '226_LO' = V6)

data_local = unite(data_local, 'ID', c('gene', 'family', 'class'), sep = ':')
row.names(data_local) = data_local$ID
data_local = select(data_local, -ID)

#################################################################
# Count table produced by TE_transcipts
#################################################################

data = read.table("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/hi_vs_lo.cntTable",header=T,row.names=1)
colnames(data) = c('214_HI', '221_HI', '226_HI', '214_LO', '221_LO', '226_LO')

TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

min_read = 1
data_transcripts = TE_data[apply(TE_data,1,function(x){max(x)}) > min_read,]

#################################################################
# Count table produced by TE_transcipts
#################################################################

differential_expression = function(count_table, p_value_cutoff = 0.05, log_fold_change_cutoff = 0.58){
  
  ## Define sampleInfo
  
  ID = colnames(count_table)
  sampleInfo = data.frame(ID,row.names=colnames(count_table))
  sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'group'), sep = '_'))
  sampleInfo$patient = factor(sampleInfo$patient)
  
  ## Construct DESeq dataset object
  
  dds <- DESeqDataSetFromMatrix(countData = count_table, colData = sampleInfo, design = ~patient + group)
  dds$group = relevel(dds$group,ref="LO")
  
  ## Run differential expression analysis
  
  dds <- DESeq(dds)
  res <- results(dds,independentFiltering=T)
  
  ## Convert results to dataframe and add signficance label
  
  results_df = as.data.frame(res)
  
  results_df = cbind(ID = rownames(results_df), results_df)
  results_df = separate(data = results_df, col = 'ID', into = c('gene', 'family', 'class'), sep = ':')
  results_df = cbind(ID = rownames(results_df), results_df)
  
  results_df = mutate(results_df, significant = case_when(padj < p_value_cutoff ~ TRUE, padj >= p_value_cutoff ~ FALSE))
  
  results_df = mutate(results_df, FC_significant = case_when(abs(log2FoldChange) > log_fold_change_cutoff ~ TRUE, 
                                                             abs(log2FoldChange) <= log_fold_change_cutoff ~ FALSE))
  
  results_df = mutate(results_df, overall_significant = case_when((significant == TRUE) & (FC_significant == TRUE) ~ TRUE, 
                                                                  (significant == FALSE) | (FC_significant == FALSE) ~ FALSE))
  
  results_df = mutate(results_df, abs_log2FoldChange = abs(log2FoldChange))
  
  sig_results_df = results_df[results_df$significant == TRUE,]
  
  ## Transform raw count data 
  
  vs_dds <- vst(dds, blind=FALSE)
  transformed_counts = as.data.frame(assay(vs_dds))
  
  sigGenes = rownames(results_df[results_df$significant == TRUE,])
  sig_transformed_counts = assay(vs_dds)[rownames(transformed_counts) %in% sigGenes,]
  
  upGenes = rownames(results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange > 0),])
  
  downGenes = rownames(results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange < 0),])
  
  return(dds)
  
}

dds = differential_expression(data_local)

x = differential_expression(data_transcripts)
combined_results_df = merge(x, y, by = 'ID')

#################################################################
# Plot
#################################################################

ggplot(data = combined_results_df, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
  geom_point() +
  geom_line(x = y)
