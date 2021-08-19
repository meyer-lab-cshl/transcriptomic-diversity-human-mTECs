library(edgeR)
library(tidyverse)

print('Importing variables...')

counts_annotated = readRDS('/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/counts_annotated')
ereMAP_loci = readRDS('/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/ereMAP_loci')

print('Building tissue list...')

counter = 0
tissue = vector()
for (name in colnames(counts_annotated)){
  
  if (grepl(pattern = 'bam', x = name) == T){
    
    counter = counter + 1
    tissue[counter] = stringr::str_split(string = name, pattern = '_')[[1]][2]
    
  }
  
}

print('Building DGEList...')

tissue = factor(tissue)
y = DGEList(counts = counts_annotated[, 6: 162], group = tissue, genes = counts_annotated[, c(1:5, 163:168)])

print('Calculating normalization factors...')

y = calcNormFactors(y)
design = model.matrix(~0+tissue)

#print('Estimating dispersion...')
#
#y = estimateDisp(y, design, robust = T)

print('Estimating RPKM values...')

RPKM_values = rpkm(y, log = F, gene.length = y$genes$width)

mean_RPKM = data.frame(row.names = row.names(RPKM_values))
for (i in unique(tissue)){
  
  mean_RPKM[i] = rowMeans(RPKM_values[, which(tissue %in% i)])
  
}

ereMAP_mean_RPKM = mean_RPKM[y$genes$Geneid %in% ereMAP_loci, ]

print('Saving output...')

saveRDS(ereMAP_mean_RPKM, file = '/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/ereMAP_mean_RPKM')

print('Finished')
