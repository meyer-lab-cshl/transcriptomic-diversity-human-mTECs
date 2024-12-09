################################################################################
# edgeR for differential expression analysis of TE_count/TE_local output. Also
# useful for obtaining RPKM values from raw counts.
################################################################################

library(edgeR)
library(tidyverse)
library(GenomicRanges)
library(fgsea)
library(glue)

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

## TE_annotation

TE_annotation = read.table(file = glue::glue("{working_directory}/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt"), 
                        header = 1)

TE_annotation = tidyr::separate(TE_annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':') %>%
  tidyr::separate(start.stop, into = c('start', 'end'), sep = '-') %>%
  dplyr::rename(locus = TE)

TE_annotation = makeGRangesFromDataFrame(TE_annotation, keep.extra.columns = T)
TE_annotation$width = GenomicRanges::width(TE_annotation) / 1000
TE_annotation = as.data.frame(TE_annotation)
TE_annotation = select(TE_annotation, -width)

## Gene annotation

gene_annotation = read.table(file = glue('{working_directory}/annotation_tables/gencode.v38_gene_annotation_table.txt'), 
                                         header = 1)
gene_annotation = dplyr::select(gene_annotation, 
                                c('Geneid', 'Chromosome', 'Start', 'End', 'Strand', 'Length')) %>%
  mutate(Length = Length / 1000) %>%
  dplyr::rename(seqnames = Chromosome, start = Start, end = End, strand = Strand, locus = Geneid, width.1 = Length)

## Bind into one annotation

annotation = bind_rows(gene_annotation, TE_annotation)

################################################################################
# TE local
################################################################################

## Imports raw counts from TE_local. Column names should be in the format
## 'uniqueID_tissue_batch.bam' e.g. pt214_mTEC.hi_new_Aligned.out.bam

count_tables = list.files(path=glue::glue('{working_directory}/count_tables/TE_local/Test'), 
                          pattern="*cntTable", 
                          full.names=TRUE, 
                          recursive=FALSE)

counter = 0
for (table in count_tables){
  
  print(table)
  
  input = read.table(file = table,
                     header = T,
                     row.names = 1)
#  input = extract_subset(input = input, mode = 'ereMAP') 
  input$Geneid = row.names(input)
  if (counter == 0){counts = input}
  else{counts = merge(counts, input, by = 'Geneid')}
  counter = counter + 1
  remove(input)
  
}

counts = separate(counts, col = Geneid, 
                  into = c('locus', 'gene', 'family', 'class'), 
                  sep = ':', 
                  remove = F)

counts_annotated =  merge(counts, annotation, by = 'locus')

## Generates a DGE object

counter = 0
tissue = vector()
for (name in colnames(counts_annotated)){
  
  if (grepl(pattern = 'bam', x = name) == T){
    
    counter = counter + 1
    tissue[counter] = stringr::str_split(string = name, pattern = '_')[[1]][2]
    
  }
  
}

tissue = factor(tissue)
y = DGEList(counts = counts_annotated[, 6:(counter + 5)], 
            group = tissue, 
            genes = counts_annotated[, c(1:5, (counter+6):(counter+10))])
y = calcNormFactors(y)

## Calculates mean RPKM

RPKM_values = rpkm(y, log = F, gene.length = y$genes$width.1)

mean_RPKM = data.frame(row.names = row.names(RPKM_values))
for (i in unique(tissue)){
  
  mean_RPKM[i] = rowMeans(RPKM_values[, which(tissue %in% i)])
  
}

mean_RPKM = cbind(mean_RPKM, y$genes)

ereMAP_loci = loadRDS('~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/ereMAP_loci')

ereMAP_loci_name = vector()
counter = 0
for (entry in ereMAP_loci){
  
  counter = counter + 1
  ereMAP_loci_name[counter] = stringr::str_split(string = entry, pattern = ':')[[1]][1]
  
}

ereMAP_mean_RPKM = mean_RPKM[rownames(mean_RPKM) %in% ereMAP_loci_name, ]

################################################################################
# TE_transcripts
################################################################################

count_tables = list.files(path=glue::glue('{working_directory}/count_tables/TE_transcripts'), 
                          pattern="*cntTable", 
                          full.names=TRUE, 
                          recursive=FALSE)

counter = 0
for (table in count_tables){
  
  print(table)
  
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

counts = separate(counts, col = Geneid, 
                  into = c('gene', 'family', 'class'), 
                  sep = ':', 
                  remove = F)

counter = 0
tissue = vector()
for (name in colnames(counts)){
  
  if (grepl(pattern = 'bam', x = name) == T){
    
    counter = counter + 1
    tissue[counter] = stringr::str_split(string = name, pattern = '_')[[1]][2]
    
  }
  
}

tissue = factor(tissue)
y = DGEList(counts = counts[, 5:(counter + 4)], 
            group = tissue)
y = calcNormFactors(y)

## CPM

CPM_values = cpm(y)
mean_CPM = data.frame(row.names = row.names(CPM_values))
for (i in unique(tissue)){
  
  mean_CPM[i] = rowMeans(CPM_values[, which(tissue %in% i)])
  
}


## Heatmap

my_heatmap = pheatmap(mean_RPKM, 
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

################################################################################
# Differential gene expression analysis
#
# Uses a paired design to compare TE expression across mTEC-HI
# and -LO samples.
################################################################################

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
