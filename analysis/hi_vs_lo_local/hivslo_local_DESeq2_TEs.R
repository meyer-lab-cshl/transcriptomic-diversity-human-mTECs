setwd('/Users/mpeacey/TE_thymus/analysis/hi_vs_lo_local')

## Import count matrix

data = read.table("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo_local/count_table/TE_local_hi_vs_lo.cntTable",
                  header=T,row.names=1)
colnames(data) = c('214_HI', '214_LO', '221_HI', '221_LO', '226_HI', '226_LO')