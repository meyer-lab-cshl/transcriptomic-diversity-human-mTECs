#################################################################
# Differential expression with DESeq2
#################################################################

## Set parameters

p_value_cutoff = 0.05
log_fold_change_cutoff = 0.58

## Import count matrix

data = read.table("hi_vs_lo.cntTable",header=T,row.names=1)
colnames(data) = c('214_HI', '221_HI', '226_HI', '214_LO', '221_LO', '226_LO')

## Subset into TE count matrix 'TE_data'

TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

## Filter count matrix to exclude non-expressed genes

min_read = 1
TE_data = TE_data[apply(TE_data,1,function(x){max(x)}) > min_read,]

## Define sampleInfo

ID = colnames(TE_data)
sampleInfo = data.frame(ID,row.names=colnames(TE_data))
sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'group'), sep = '_'))

## Construct DESeq dataset object

dds <- DESeqDataSetFromMatrix(countData = TE_data, colData = sampleInfo, design = ~patient + group)
dds$group = relevel(dds$group,ref="LO")

## Run differential expression analysis

dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)

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