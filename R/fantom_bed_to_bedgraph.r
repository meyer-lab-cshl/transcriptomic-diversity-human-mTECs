library(plyr)

setwd("/home/stroemic/hiwi_16/data/external/Fantom5/CTSS_BED/tissues")

fantoms <- list.files(pattern="*.bed.gz")
for(i in fantoms){
  dat <- read.table(i)
  #get identifier to name later file
  name= paste((head(unlist(strsplit(i, "[.]" )), n=-2)), collapse='.')
  print(name)
  #rename columns
  dat<- rename(dat,c("V1"="chromosome","V2"="start","V3"="end","V4"="name","V5"="count","V6"="strand"))
  #normalise counts based on total read count
  dat$count = dat$count/(sum(dat$count))*10000000
  #count negative when on negative strand
  dat[dat$strand == "-",]$count <- (dat[dat$strand == "-",]$count)*-1

  ###export dataframe from dat
  dat_export <- data.frame("chrom"=dat$chromosome, "start"=dat$start, "end"=dat$end, "dataValue"=dat$count)
  dat_export <- dat_export[complete.cases(dat_export),]
  dat_export <- dat_export[with(dat_export, order(chrom, start)),]
  #export data frame as csv, do not use quotes, rownames; do use colnames
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/data/Summary_counts/bedgraph/",name, ".bedgraph", sep=""), row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
  print(paste("/home/stroemic/hiwi_16/data/Summary_counts/bedgraph/",name, ".bedgraph", sep=""))
  
}
