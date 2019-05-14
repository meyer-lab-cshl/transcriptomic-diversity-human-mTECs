library(plyr)

setwd("/home/stroemic/hiwi_16/analysis/STAR_alignment_hg19/5Seq_processed_IVT/")

files <- grep("collapsing_1",list.files(pattern = "_countsPerMol.tsv$", recursive = TRUE), value=TRUE)
for(i in files){
  dat <- read.table(i, header=TRUE)
  #get identifier to name later file
  name= paste(head(unlist(strsplit(unlist(strsplit(i, "/"))[3], "_")), n=-1), collapse='_')
  print(name)
  #drop IVT
  dat <- dat[!(dat$chr == "pGIBS-LYS" | dat$chr == "pGIBS-PHE" | dat$chr == "pGIBS-THR"),]
  #drop molbc and calculate real count based on positions
  dat <- subset(dat, select = -c(molbc))
  dat$count = 1
  dat<- aggregate(dat$count, by=list(dat$chr, dat$pos, dat$strand), FUN=sum)
  #rename columns
  dat<- rename(dat,c("Group.1"="chromosome","Group.2"="start","Group.3"="strand", "x"="count"))
  dat$end <- as.integer(dat$start +1)
  #normalise counts based on total read count
  dat$count = dat$count/(sum(dat$count))*10000000
  dat[dat$strand == "-",]$count <- (dat[dat$strand == "-",]$count)*-1
  dat <- dat[with(dat, order(chromosome, start)),]
  ###export dataframe from dat
  dat_export <- data.frame("chrom"=dat$chromosome, "start"=dat$start, "end"=dat$end, "dataValue"=dat$count)
  dat_export <- dat_export[complete.cases(dat_export),]
  #export data frame as csv, do not use quotes, rownames; do use colnames
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/data/Summary_counts/bedgraph/",name, ".bedgraph", sep=""), row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
  print(paste("/home/stroemic/hiwi_16/data/Summary_counts/bedgraph/",name, ".bedgraph", sep=""))
}