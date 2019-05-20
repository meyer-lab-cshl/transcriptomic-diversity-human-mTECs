library(plyr)
library(data.table)

setwd("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/")

files <- list.files(pattern="*Muenster.filtered.csv")

datalist=list()
c <- 1
for(i in files){
  print(i)
  dat <- read.table(i, header = TRUE, sep = ",")
  dat$iteration <- c
  datalist[[c]] <- dat
  c <- c+1
}

dat_internal = do.call(rbind, datalist)

#use either sum or mean as function
dat_agg_internal <- setDT(dat_internal)[, lapply(.SD, mean), by=.(chrom, geneID, strand, start, end ), .SDcols=c("BarcodeCount")]
setDF(dat_agg_internal)

write.table(dat_agg_internal, file ="/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered-mean.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
