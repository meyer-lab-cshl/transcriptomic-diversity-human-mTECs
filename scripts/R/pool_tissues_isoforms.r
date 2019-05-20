library(plyr)
library(data.table)

setwd("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/")

files <- list.files(pattern="*ctss.isoforms.csv")

datalist<-list()
datalist_wo<-list()
c <- 1
for(i in files){
  print(i)
  dat <- read.table(i, header = TRUE, sep = ",")
  dat$iteration <- c
  if (grepl("thymus", i)) {
    print("one")
    datalist[[c]] <- dat
  } else {
    print("both")
    datalist[[c]] <- dat
    datalist_wo[[c]] <- dat
  }
  c <- c+1
}
print("rbinding")
dat_fantom5 = do.call(rbind, datalist)
dat_fantom5_wo_thymus = do.call(rbind, datalist_wo)

print("aggregating")
##either do lapply(.SD, sum) or (.SD, mean)
dat_agg_fantom5 <- setDT(dat_fantom5)[, lapply(.SD, mean), by=.(geneID, Position), .SDcols=c("BarcodeCount","ReadCount","PosFromAnno","Class")]
setDF(dat_agg_fantom5)
dat_agg_fantom5_wo_thymus <- setDT(dat_fantom5_wo_thymus)[, lapply(.SD, mean), by=.(geneID, Position), .SDcols=c("BarcodeCount","ReadCount","PosFromAnno","Class")]
setDF(dat_agg_fantom5_wo_thymus)

print("writing out")
##change file extension sum or mean 
write.table(dat_agg_fantom5, file ="/home/stroemic/hiwi_16/analysis/lars_perl/all_63/all.tissues.isoforms_mean.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(dat_agg_fantom5_wo_thymus, file ="/home/stroemic/hiwi_16/analysis/lars_perl/all_63/all.tissues.wo.thymus.isoforms_mean.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

