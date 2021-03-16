library(GenomicRanges)
library(dplyr)
library(tidyr)

TE_annotation = read.table(file = "/Users/mpeacey/TE_thymus/analysis/hg38_rmsk_TE.gtf.locInd.locations.txt", header = 1)
TE_annotation = separate(TE_annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':')
TE_annotation = separate(TE_annotation, start.stop, into = c('start', 'end'), sep = '-')
TE_annotation = rename(TE_annotation, locus = TE)

df = merge(results_df_local, TE_annotation, by = 'locus')

GRanges_TE = makeGRangesFromDataFrame(df, keep.extra.columns = T)