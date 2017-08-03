library(plyr)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

setwd("/home/stroemic/hiwi_16/analysis/lars_perl/mouse/")

my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chrM", "chrX", "chrY")

#####Load gene object
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)

isoforms <- list.files(pattern=".*-.*.isoforms.csv")
start.time <- Sys.time()
for(i in isoforms){
  dat <- read.table(i, header=TRUE, sep = ",")
  #get identifier to name later file
  if (grepl("_", i)) {
    name= paste((head(unlist(strsplit(i, "[_]" )), n=-1)), collapse='_')
    print(name)
  } else {
    name= paste((head(unlist(strsplit(i, "[.]" )), n=-2)), collapse='.')
    print(name)
  }
  #rename columns
  dat<- rename(dat,c("Position"="start","BarcodeCount"="count"))
  dat <- subset(dat, select = -c(ReadCount, PosFromAnno, Class))
  ###split data frame in genes and intergenic regions
  dat_inter <- dat[grepl("intergenic", dat$geneID),]
  dat <- dat[!grepl("intergenic", dat$geneID),]
  ###for genes annotate chromosome and strand with tx_genes object
  dat$chrom <- with(dat, as.character(seqnames(tx_genes[as.character(dat$geneID)])) )
  dat$strand <- with(dat, as.character(strand(tx_genes[dat$geneID])) )
  ###for intergenic regions annotate chromosome from name and strand as plus
  dat_inter$chrom <- gsub(".*_(chr.{1,2})_.*", "\\1", dat_inter$geneID )
  dat_inter$strand <- gsub(".*_.*_(.)_.*", "\\1", dat_inter$geneID )
  ###merge dataframes again
  dat <- rbind(dat, dat_inter)
  ###define end, adjust count based on strand, resort dataframe
  dat$end <- as.integer(dat$start +1)
  dat[dat$strand == "-",]$count <- (dat[dat$strand == "-",]$count)*-1
  dat <- dat[with(dat, order(chrom, start)),]
 
  ###export dataframe from dat
  dat_export <- data.frame("chrom"=dat$chrom, "start"=dat$start, "end"=dat$end, "dataValue"=dat$count)
  dat_export <- dat_export[complete.cases(dat_export),]
  dat_export <- dat_export[with(dat_export, order(chrom, start)),]
  #export data frame as csv, do not use quotes, rownames; do use colnames
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/analysis/lars_perl/mouse/",name, ".isoforms.bedgraph", sep=""), row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
  print(paste("/home/stroemic/hiwi_16/analysis/lars_perl/mouse/",name, ".isoforms.bedgraph", sep=""))
  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


