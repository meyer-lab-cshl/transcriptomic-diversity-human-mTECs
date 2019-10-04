library(plyr)
library(GenomicRanges)

setwd("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/")

files <- list.files(pattern = "\\.filtered.csv$")
for(i in files){
  #####consensus sites csv read in
  dat <- read.table(i, sep = ",", header = TRUE)
  #get identifier to name later file
  name= paste(head(unlist(strsplit(i, "[.]")), n=-2), collapse='.')
  #name= paste(head(unlist(strsplit(i, "_")), n=-1), collapse='_')
  print(name)
    dat[dat$strand == "-",]$BarcodeCount <- (dat[dat$strand == "-",]$BarcodeCount)*-1
  dat <- dat[with(dat, order(chrom, start)),]
  
  dat_export <- data.frame("chrom"=dat$chrom, "start"=dat$start, "end"=dat$end, "dataValue"=dat$BarcodeCount)
  dat_export <- dat_export[complete.cases(dat_export),]
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/bedgraph/",name, ".filtered.bedgraph", sep=""), row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
}
