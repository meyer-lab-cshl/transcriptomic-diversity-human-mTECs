#################################################################
# GTEX_test
#################################################################

## Make count table

data = read.table("/Users/mpeacey/TE_thymus/analysis/count_tables/TE_local_GTEX_test",
                  header=T,row.names=1)
colnames(data) = c('skin_1', 'skin_2', 'brain_1', 'brain_2')

min_read = 1
data = data[apply(data,1,function(x){max(x)}) > min_read,]

# Run differential expression

ID = colnames(data)
sampleInfo = data.frame(ID,row.names=colnames(data))
sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('tissue', 'patient'), sep = '_'))
sampleInfo$patient = factor(sampleInfo$patient)

## Construct DESeq dataset object

dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~tissue)
dds$tissue = relevel(dds$tissue,ref="skin")

## Run differential expression analysis

dds = DESeq(dds)
res = results(dds, independentFiltering = F)

results_df = as.data.frame(res)

results_df_gene = results_df[grep("^ENSG",rownames(results_df)),]