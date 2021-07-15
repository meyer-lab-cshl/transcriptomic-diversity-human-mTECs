library(DESeq2)
library(ggplot2)
library(tidyverse)
library(glue)

functions_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_functions/"
functions = c('extract_subset', 'differential_expression')

for (i in functions){
  
  load(glue::glue('{functions_directory}{i}'))
  
}

#################################################################
# DESeq2
#################################################################

## 'count_table_directory' should contain the text file 'TE_transcripts_counts' containing
## the raw counts from each tissue of interest. Each entry should be labelled in the format 
## '{unique ID}_{tissue}_{batch}'. e.g. 'pt214_mTEC-hi_our-data'

## Data import:

count_table_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/"
transcripts_data = read.table(glue::glue('{count_table_directory}TE_transcripts_counts'),header=T,row.names=1)

## Run DESeq2 to obtain normalized counts

TE_data = extract_subset(mode = 'TE', input = transcripts_data)

dds_transcripts_TE = differential_expression(TE_data, design=~tissue)

vs_dds_transcripts_TE = vst(dds_transcripts_TE, blind=FALSE)

assay(vs_dds_transcripts_TE) = limma::removeBatchEffect(assay(vs_dds_transcripts_TE), vs_dds_transcripts_TE$batch)

#################################################################
# PCA w/ GTEx data (A)
#################################################################

## Plot

pcaData = plotPCA(vs_dds_transcripts_TE, intgroup='tissue', returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

## Use this to have only a subset of tissues in color: useful when you have lots of tissues
## and you're only interested in a few.

colored_tissues = c('mTEC.hi', 'mTEC.lo', 'ESC', 'Testis')

pcaData = mutate(pcaData, color = case_when(tissue %in% colored_tissues ~ T,
                                   !(tissue %in% colored_tissues) ~ F))

PCA = ggplot(pcaData, aes(PC1, PC2)) + 
  geom_point(data = pcaData, aes(x = PC1, y = PC2), color = '#9A9A9A', size = 3) +
  geom_point(data = subset(pcaData, color == F), size=3, shape = 21, stroke = 0) +
  geom_point(data = subset(pcaData, color == T), aes(fill = tissue), size=3, shape = 21, stroke = 0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  labs(fill= "Tissue") +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff', '#55A257', '#E93C00'))

PCA + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                         plot.subtitle = element_text(size = 14),
                         axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),
                         axis.title = element_text(size = 14),
                         axis.line = element_line(size = 0.8),
                         panel.border = element_blank(),
                         legend.text = element_text(size = 15),
                         legend.title = element_text(size = 0),
                         legend.position="top",
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/A_pca.png", 
       width = 5.25, height = 6, units = "in")

#################################################################
# Fraction of reads mapping to TEs (supplement A?)
#################################################################

dds_transcripts = differential_expression(transcripts_data, design=~tissue)

df_gene_counts = data.frame(gene_counts = colSums(counts(extract_subset(mode = 'gene', input = dds_transcripts), normalized = T)))
df_gene_counts = cbind(ID = rownames(df_gene_counts), df_gene_counts)

df_TE_counts = data.frame(TE_counts = colSums(counts(extract_subset(mode = 'TE', input = dds_transcripts), normalized = T)))
df_TE_counts = cbind(ID = rownames(df_TE_counts), df_TE_counts)

df = merge(df_gene_counts, df_TE_counts, by = 'ID')

df = mutate(df, total = gene_counts + TE_counts)
df = mutate(df, percent_TEs = (TE_counts / total)*100)

df = separate(df, col = ID, into = c('patient', 'tissue'), sep = '_')

df = mutate(df, tissue = forcats::fct_reorder(tissue, percent_TEs))

## Boxplot

df = mutate(df, color = case_when(tissue %in% colored_tissues ~ T,
                                      !(tissue %in% colored_tissues) ~ F))

plot = ggplot(data = df, aes(x = tissue, y = percent_TEs)) +
  geom_boxplot(data = subset(df, color == F), fill = '#9A9A9A', width = 0.6) +
  geom_boxplot(data = subset(df, color == T), aes(fill = tissue), width = 0.6) +
  geom_point() +
  ylab('% of reads mapping to TEs') +
  xlab('Tissue/cell') +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff', '#55A257', '#E93C00'))

plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 16),
                          plot.subtitle = element_text(size = 12),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_line(color = 'gray'),
                          axis.text.x = element_blank(),
                          axis.text.y = element_text(size = 13),
                          axis.title.x = element_text(size = 15, margin = margin(t = 12.5)),
                          axis.title.y = element_text(size = 15, margin = margin(r = 7.5)),
                          axis.line = element_line(size = 0.8),
                          panel.border = element_blank(),
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size = 14))

ggsave("~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/overall_TE_abundance.png", 
       width = 7, height = 5, units = "in")
