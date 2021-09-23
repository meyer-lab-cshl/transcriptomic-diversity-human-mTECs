################################################################################
# Get RPKM counts from TE_local output
################################################################################

library(edgeR)
library(tidyverse)
library(glue)

working_directory = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis'
functions_directory = glue("{working_directory}/R_functions/")

## Function import
functions = c('extract_subset')

for (i in functions){
  
  load(glue('{functions_directory}{i}'))
  
}

## Variable import
ereMAP_loci = readRDS('{working_directory}/R_variables/ereMAP_loci')

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
  dplyr::rename(locus = TE)

annotation = makeGRangesFromDataFrame(annotation, keep.extra.columns = T)
annotation$width = GenomicRanges::width(annotation)
annotation = as.data.frame(annotation)

################################################################################
# Prepare raw counts
#
# Imports raw counts from TE_local. Column names should be in the format
# 'uniqueID_tissue_batch.bam' e.g. pt214_mTEC.hi_new_Aligned.out.bam
################################################################################

count_tables = list.files(path=glue::glue('{working_directory}/count_tables/TE_local/New'), 
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
                  into = c('locus', 'gene', 'family', 'class'), 
                  sep = ':', 
                  remove = F)

counts_annotated =  merge(counts, annotation, by = 'locus')

################################################################################
# Get RPKM values
################################################################################

counter = 0
tissue = vector()
for (name in colnames(counts_annotated)){
  
  if (grepl(pattern = 'bam', x = name) == T){
    
    counter = counter + 1
    tissue[counter] = stringr::str_split(string = name, pattern = '_')[[1]][2]
    
  }
  
}

tissue = factor(tissue)
y = DGEList(counts = counts_annotated[, 6: 162], group = tissue, genes = counts_annotated[, c(1:5, 163:168)])
y = calcNormFactors(y)
design = model.matrix(~0+tissue)

## Step takes a long time and not sure if necessary
#y = estimateDisp(y, design, robust = T)

RPKM_values = rpkm(y, log = F, gene.length = y$genes$width)

# Mean RPKM + heatmap

mean_RPKM = data.frame(row.names = row.names(RPKM_values))
for (i in unique(tissue)){
  
  mean_RPKM[i] = rowMeans(RPKM_values[, which(tissue %in% i)])
  
}

ereMAP_mean_RPKM = mean_RPKM[y$genes$Geneid %in% ereMAP_loci, ]

saveRDS(ereMAP_mean_RPKM, file = glue('{working_directory}/R_variables/ereMAP_mean_RPKM'))

#my_heatmap = pheatmap(ereMAP_mean_RPKM, 
#                      cluster_rows=T,
#                      show_rownames=F,
#                      show_colnames = T,
#                      cluster_cols=T,
#                      scale = 'row',
#                      angle_col = 45)

