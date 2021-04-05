library(DESeq2)
library(dplyr)
library(tidyr)
library(glue)
library(stringr)
library(limma)

#################################################################
# Functions
#################################################################

## Differential expression

differential_expression = function(raw_count_table, min_reads = 2, design = ~tissue){
  
  ## Define sampleInfo
  
  ID = colnames(raw_count_table)
  sampleInfo = data.frame(ID,row.names=colnames(raw_count_table))
  sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'tissue', 'batch'), sep = '_'))
  sampleInfo$patient = factor(sampleInfo$patient)
  
  ## Construct DESeq dataset object
  
  dds <- DESeqDataSetFromMatrix(countData = raw_count_table, colData = sampleInfo, design = design)
  #dds$group = relevel(dds$group,ref="LO")
  
  ## Run differential expression analysis
  
  dds = DESeq(dds)
  
  ## Pre-filter to remove rows with less than 'min_reads' total (default 2)
  
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

standardize_column_names = function(raw_counts){
  
  new_col_name = vector(, length(colnames(raw_counts)))
  
  for (i in 1:length(colnames(raw_counts))){
    
    if (colnames(raw_counts)[i] %in% colnames(mTEC_counts)){
      
      a = str_split(colnames(raw_counts)[i], '_')[[1]][1]
      
      b = 'mTEC'
      c = str_split(colnames(raw_counts)[i], '_')[[1]][2]
      b = paste(b, c, sep = '-')
      b = paste(b, 'new', sep = '_')
      
    }
    
    else{
      
      a = str_split(colnames(raw_counts)[i], '_')[[1]][1]
      
      if (colnames(raw_counts)[i] %in% colnames(testis_counts)){
        
        b = 'Testis_GTEx'
        
      }
      
      if (colnames(raw_counts)[i] %in% colnames(ovaries_counts)){
        
        b = 'Ovary_GTEx'
        
      }
      
      if (colnames(raw_counts)[i] %in% colnames(muscle_counts)){
        
        b = 'Muscle_GTEx'
        
      }
      
      if (colnames(raw_counts)[i] %in% colnames(ESC_counts)){
        
        b = 'ESC_UCSC'
        
      }
      
    }
    
    new_col_name[i] = paste(a, b, sep = '_')
    
  }
  
  colnames(raw_counts) = new_col_name
  
  return(raw_counts)
  
}

#################################################################
# TE transcripts
#################################################################

## Data import

raw_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_count/TE_transcripts_counts",header=T,row.names=1)
testis_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_count/testis_counts",header=T,row.names=1)
ovaries_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_count/ovaries_counts",header=T,row.names=1)
mTEC_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_count/TE_transcripts_hi_vs_lo.cntTable",header=T,row.names=1)
muscle_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_count/muscle_counts",header=T,row.names=1)
ESC_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_count/ESC_counts",header=T,row.names=1)

data = standardize_column_names(raw_counts = raw_counts)

## Run DESeq2

dds_transcripts = differential_expression(data, design=~tissue)
dds_transcripts_gene = extract_from_DESeq2(mode = 'gene', input = dds_transcripts)
dds_transcripts_TE = extract_from_DESeq2(mode = 'TE', input = dds_transcripts)

## Normalized counts

vs_dds_transcripts_gene = vst(dds_transcripts_gene, blind=FALSE)
assay(vs_dds_transcripts_gene) = limma::removeBatchEffect(assay(vs_dds_transcripts_gene), vs_dds_transcripts_gene$batch)

vs_dds_transcripts_TE = vst(dds_transcripts_TE, blind=FALSE)
assay(vs_dds_transcripts_TE) = limma::removeBatchEffect(assay(vs_dds_transcripts_TE), vs_dds_transcripts_TE$batch)

## Differential expression

results_transcripts = results(dds_transcripts, 
                              contrast = c('tissue', 'mTEC-hi', 'mTEC-lo'), 
                              independentFiltering = F)

results_transcripts_gene = extract_from_DESeq2(mode = 'gene', input = results_transcripts)
results_transcripts_TE = extract_from_DESeq2(mode = 'TE', input = results_transcripts)

results_df_transcripts_gene = process_DESeq2_results(results = results_transcripts_gene, mode = 'Gene')
results_df_transcripts_TE = process_DESeq2_results(results = results_transcripts_TE, mode = 'TE_transcripts')

results_df_transcripts_gene_sigdiff = filter(results_df_transcripts_gene, significant == T)

results_df_transcripts_TE_sigdiff = filter(results_df_transcripts_TE, significant == T)

#################################################################
# TE_local
#################################################################

## Data import

raw_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_local/TE_local_counts",header=T,row.names=1)

mTEC_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_local/TE_local_hi_vs_lo.cntTable_old",header=T,row.names=1)
testis_counts = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_local/testis_local_count_table",header=T,row.names=1)

data = standardize_column_names(raw_counts = mTEC_counts)

## Run DESeq2

dds_local = differential_expression(data, design=~patient+tissue)
dds_local_gene = extract_from_DESeq2(mode = 'gene', input = dds_local)
dds_local_TE = extract_from_DESeq2(mode = 'TE', input = dds_local)

## Normalized counts

vs_dds_local_gene = vst(dds_local_gene, blind=FALSE)
assay(vs_dds_local_gene) = limma::removeBatchEffect(assay(vs_dds_local_gene), vs_dds_local_gene$batch)
#assay(vs_dds_local_gene) = limma::removeBatchEffect(assay(vs_dds_local_gene), vs_dds_local_gene$patient)

vs_dds_local_TE = vst(dds_local_TE, blind=FALSE)
assay(vs_dds_local_TE) = limma::removeBatchEffect(assay(vs_dds_local_TE), vs_dds_local_TE$batch)
#assay(vs_dds_local_TE) = limma::removeBatchEffect(assay(vs_dds_local_TE), vs_dds_local_TE$patient)

## Differential expression

results_local = results(dds_local, 
                        independentFiltering = F)

results_local = results(dds_local, 
                        independentFiltering = F,
                        contrast = c('tissue', 'mTEC-hi', 'mTEC-lo'))

results_local_gene = extract_from_DESeq2(mode = 'gene', input = results_local)
results_local_TE = extract_from_DESeq2(mode = 'TE', input = results_local)

results_df_local_gene = process_DESeq2_results(results = results_local_gene, mode = 'Gene')
results_df_local_TE = process_DESeq2_results(results = results_local_TE, mode = 'TE_local')

results_df_local_gene_up = filter(results_df_local_gene, (significant == T) & (log2FoldChange > 0))
results_df_local_gene_unchanged = filter(results_df_local_gene, significant == F)
results_df_local_gene_down = filter(results_df_local_gene, (significant == T) & (log2FoldChange < 0))
results_df_local_gene_sigdiff = filter(results_df_local_gene, significant == T)

results_df_local_TE_up = filter(results_df_local_TE, (significant == T) & (log2FoldChange > 0))
results_df_local_TE_unchanged = filter(results_df_local_TE, significant == F)
results_df_local_TE_down = filter(results_df_local_TE, (significant == T) & (log2FoldChange < 0))
results_df_local_TE_sigdiff = filter(results_df_local_TE, significant == T)

#################################################################
# TE_local (without locus information)
#################################################################

data = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_local_hi_vs_lo.cntTable",
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
