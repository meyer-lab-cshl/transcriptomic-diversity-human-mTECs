library(plyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

setwd("/home/stroemic/hiwi_16/analysis/lars_perl/all_63")
my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

#####Load gene object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)

files <- list.files(pattern = "isoforms.csv$")
for(i in files){
  #####consensus sites csv read in
  dat <- read.table(files[1], sep = ",", header = TRUE)
  #get identifier to name later file
  name= paste(head(unlist(strsplit(i, "[.]")), n=-2), collapse='.')
  #name= paste(head(unlist(strsplit(i, "_")), n=-1), collapse='_')
  print(name)
  
  dat <- subset(dat, select = -c(ReadCount, PosFromAnno, Class))
  dat <- dat[!grepl("intergenic", dat$geneID),]
  dat$chrom <- with(dat, as.character(seqnames(tx_genes[dat$geneID])) )
  dat$strand <- with(dat, as.character(strand(tx_genes[dat$geneID])) )
  dat$end <- as.integer(dat$Position +1)
  dat[dat$strand == "-",]$BarcodeCount <- (dat[dat$strand == "-",]$BarcodeCount)*-1
  dat <- dat[with(dat, order(chrom, Position)),]
  
  dat_export <- data.frame("chrom"=dat$chrom, "start"=dat$Position, "end"=dat$end, "dataValue"=dat$BarcodeCount)
  dat_export <- dat_export[complete.cases(dat_export),]
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/bedgraph/",name, ".isoforms.bedgraph", sep=""), row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
}
