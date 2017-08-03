library(plyr)
library(data.table)

setwd("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/")

files <- list.files(pattern="*collapsed.isoforms.csv")

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
dat_agg_internal <- setDT(dat_internal)[, lapply(.SD, mean), by=.(geneID, Position), .SDcols=c("BarcodeCount","ReadCount","PosFromAnno","Class")]
setDF(dat_agg_internal)

write.table(dat_agg_internal, file ="/home/stroemic/hiwi_16/analysis/lars_perl/all_63/all_mTECs.isoforms_mean.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
