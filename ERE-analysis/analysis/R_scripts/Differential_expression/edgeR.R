################################################################################
# edgeR for differential expression analysis of TE_count/TE_local output. Also
# useful for obtaining RPKM values from raw counts.
################################################################################

library(edgeR)
library(tidyverse)
library(GenomicRanges)
library(fgsea)

working_directory = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis'
functions_directory = glue("{working_directory}/R_functions/")

functions = c('extract_subset')

for (i in functions){
  
  load(glue('{functions_directory}{i}'))
  
}

################################################################################
# Prepare annotation
# 
# Imports the RepeatMasker annotation (downloaded from 
# http://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/prebuilt_indices/)
################################################################################

annotation = read.table(file = glue::glue("{working_directory}/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt"), 
                        header = 1)

annotation = tidyr::separate(annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':') %>%
  tidyr::separate(start.stop, into = c('start', 'end'), sep = '-') %>%
  dplyr::rename(annotation, locus = TE)

annotation = makeGRangesFromDataFrame(annotation, keep.extra.columns = T)
annotation$width = GenomicRanges::width(annotation)
annotation = as.data.frame(annotation)

#################################################################
# Prepare raw counts
#
# Imports raw counts from TE_local.
#################################################################

count_tables = list.files(path='/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/TE_local/New', 
                          pattern="*cntTable", 
                          full.names=TRUE, 
                          recursive=FALSE)

counter = 0
for (table in count_tables){
  
  input = read.table(file = table,
                     header = T,
                     row.names = 1)
  input = extract_subset(input = input, mode = 'ERE') 
  input$Geneid = row.names(input)
  if (counter == 0){counts = input}
  else{counts = merge(counts, input, by = 'Geneid')}
  counter = counter + 1
  remove(input)
  
}

## mTECs

mTEC_counts = read.table('/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/TE_local/TE_local_hi_vs_lo.cntTable'
                  ,header=T,
                  row.names=1)
mTEC_ERE_counts = extract_subset(input = mTEC_counts, mode = 'ERE') 
mTEC_ERE_counts$Geneid = row.names(mTEC_ERE_counts)
mTEC_ERE_counts = separate(mTEC_ERE_counts, col = Geneid, into = c('locus', 'gene', 'family', 'class'), sep = ':', remove = F)

## ESCs

ESC_counts = read.table('/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/TE_local/ESC_TE_local.cntTable'
                        ,header=T,
                        row.names=1)
ESC_ERE_counts = extract_subset(input = ESC_counts, mode = 'ERE') 
ESC_ERE_counts$Geneid = row.names(ESC_ERE_counts)
ESC_ERE_counts = separate(ESC_ERE_counts, col = Geneid, into = c('locus', 'gene', 'family', 'class'), sep = ':', remove = F)

## Testis

testis_counts = read.table('/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/TE_local/Testis_TE_local.cntTable'
                        ,header=T,
                        row.names=1)
testis_ERE_counts = extract_subset(input = testis_counts, mode = 'ERE') 
testis_ERE_counts$Geneid = row.names(testis_ERE_counts)
testis_ERE_counts = separate(testis_ERE_counts, col = Geneid, into = c('locus', 'gene', 'family', 'class'), sep = ':', remove = F)

## All counts

all_ERE_counts = merge(mTEC_ERE_counts, ESC_ERE_counts) %>%
  merge(testis_ERE_counts)

#################################################################
# Differential gene expression analysis
#
# Uses a paired design to compare TE expression across mTEC-HI
# and -LO samples.
#################################################################

mTEC_ERE_counts_annotated =  merge(mTEC_ERE_counts, annotation, by = 'locus')

patient = factor(c('pt226', 'pt226', 'pt221', 'pt221', 'pt214', 'pt214'))
tissue = factor(c('lo', 'hi', 'lo', 'hi', 'hi', 'lo'))
y = DGEList(counts = mTEC_ERE_counts_annotated[, 2:7], group = tissue, 
            genes = mTEC_ERE_counts_annotated[, c(1, 8:15)])
y = y[filterByExpr(y), , keep.lib.sizes=FALSE]
y = calcNormFactors(y)
design = model.matrix(~0+patient+tissue)
y = estimateDisp(y, design, robust = T)

## Calculate FDR

fit = glmFit(y, design)
lrt = glmLRT(fit)

edgeR_results = topTags(lrt, n = nrow(lrt$table))$table

## Annotation

edgeR_results$ID = row.names(edgeR_results)

edgeR_results = mutate(edgeR_results, significant = case_when(FDR < 0.05 ~ T,
                                                              FDR > 0.05 ~ F))

edgeR_results = mutate(edgeR_results, ereMAP = case_when(Geneid %in% ereMAP_loci ~ T,
                                                         !(Geneid %in% ereMAP_loci) ~ F))

## Export

saveRDS(object = edgeR_results,
        file = glue::glue('{working_directory}/R_variables/edgeR_results_local_ERE'))

#################################################################
# Get RPKM values
#################################################################

ERE_counts_annotated =  merge(all_ERE_counts, annotation, by = 'locus')

tissue = factor(c('lo', 'hi', 'lo', 'hi', 'hi', 'lo', 'ESC', 'ESC', 'testis', 'testis'))
y = DGEList(counts = ERE_counts_annotated[, 6: 15], group = tissue, genes = ERE_counts_annotated[, c(1:5, 16:19)])
y = calcNormFactors(y)
design = model.matrix(~0+tissue)
y = estimateDisp(y, design, robust = T)

RPKM_values = rpkm(y, log = F, gene.length = y$genes$width)

# Mean RPKM + heatmap

mean_RPKM = data.frame(row.names = row.names(RPKM_values))
for (i in unique(tissue)){
  
  mean_RPKM[i] = rowMeans(RPKM_values[, which(tissue %in% i)])
  
}

ereMAP_mean_RPKM = mean_RPKM[y$genes$Geneid %in% ereMAP_loci, ]

my_heatmap = pheatmap(ereMAP_mean_RPKM, 
                      cluster_rows=T,
                      show_rownames=F,
                      show_colnames = T,
                      cluster_cols=T,
                      scale = 'row',
                      angle_col = 45)

# Ranked by expression in one tissue

RPKM_df = y$genes

for (i in unique(tissue)){
  
  RPKM_df[i] = mean_RPKM[, i]
  
}

RPKM_df = mutate(RPKM_df, ereMAP = case_when(Geneid %in% ereMAP_loci ~ T,
                                             !(Geneid %in% ereMAP_loci) ~ F))

RPKM_df_hi = select(RPKM_df, c('locus', 'hi', 'ereMAP')) %>%
  mutate(ID = forcats::fct_reorder(locus, hi))
RPKM_df_lo = select(RPKM_df, c('locus', 'lo', 'ereMAP'))
mutate(ID = forcats::fct_reorder(locus, lo))

filtered_RPKM = subset(RPKM_df, mean_RPKM > 0.5)

plot = ggplot() +
  geom_point(data = filtered_RPKM, aes(x = ID, y = log2(mean_RPKM)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(filtered_RPKM, ereMAP == T), aes(x = ID, y = log2(mean_RPKM), fill = ereMAP), size = 2, alpha = 0.8, shape = 21, stroke = 0) +
  geom_point(data = subset(filtered_RPKM, ereMAP == F), aes(x = ID, y = log2(mean_RPKM)), size = 1, alpha = 0.8, shape = 21, stroke = 0) +
  scale_fill_manual(values = c('#e41a1c'))

plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                          plot.subtitle = element_text(size = 14),
                          axis.text.x = element_blank(),
                          axis.text.y = element_text(size = 14),
                          axis.title = element_text(size = 14),
                          axis.line = element_line(size = 0.8),
                          panel.border = element_blank(),
                          legend.text = element_text(size = 15),
                          legend.title = element_text(size = 18),
                          legend.position = c(0.2, 0.93),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())

# GSEA

ranks = RPKM_df$mean_RPKM
names(ranks) = RPKM_df$ID
ranks = sort(ranks)

pathways = list('ereMAPs' = subset(RPKM_df, ereMAP == T)$ID)

fgseaRes = fgsea(pathways = pathways, stats = ranks)
plotEnrichment(pathways$ereMAPs, ranks)


