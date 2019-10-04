library(plyr)
library(data.table)

setwd("/home/stroemic/hiwi_16/data/Summary_counts/bed")

files <- list.files(pattern="*.bed")

datalist=list()
c <- 1
for(i in files){
  print(i)
  dat <- read.table(i)
  dat <- rename(dat,c("V1"="chrom","V2"="start","V3"="end","V4"="name", "V5"="count","V6"="strand"))
  dat$iteration <- c
  datalist[[c]] <- dat
  c <- c+1
}

dat_internal = do.call(rbind, datalist)

dat_agg_internal <- setDT(dat_internal)[, lapply(.SD, sum), by=.(chrom, start, end, name, strand), .SDcols=c("count", "iteration")]
setDF(dat_agg_internal)

write.table(dat_agg_internal[,c(1:4,6,5)], file ="/home/stroemic/hiwi_16/data/Summary_counts/bed/all_mTECs.bed", row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
