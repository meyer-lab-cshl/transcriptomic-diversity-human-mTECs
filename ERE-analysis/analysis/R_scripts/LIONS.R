library("ggvenn")
library(tidyverse)
library(GenomicRanges)

LIONS = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/External_packages/LIONS/mTEC-analysis/mTEC-analysis.rslions', 
         sep = '\t', header = T) %>%
  separate(coordinates, into = c('chromosome', 'range'), sep = ':') %>%
  separate(range, into = c('start', 'end'), sep = '-') %>%
  separate(repeatName, remove = F, into = c('sub-family', 'class', 'family'), sep = ':')

GRanges_LIONS = makeGRangesFromDataFrame(LIONS, keep.extra.columns = T)

LIONS_overlaps = findOverlaps(query = GRanges_LIONS, subject = GRanges_ERE)
LIONS_overlaps = as.data.frame(LIONS_overlaps)

for (i in 1:nrow(LIONS)){
  
  ERE = LIONS[i, 'sub-family']
  print(ERE)
  
}

venn_diagram_list = list('mTEC-LO' = subset(LIONS, Normal_Occ >= 2)$transcriptID,
                         'mTEC-HI' = subset(LIONS, Cancer_Occ >= 2)$transcriptID)

ggvenn(venn_diagram_list)

LIONS_all = subset(LIONS, Normal_Occ >= 2 | Cancer_Occ >= 2)$repeatName