library(plyr)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

setwd("/home/stroemic/hiwi_16/analysis/STAR_alignment_mm9/5Seq_processed/")

my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chrM", "chrX", "chrY")

#####Load gene object
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
tx_genes <- tx_genes + 100
o <- order(as.numeric(tx_genes$gene_id))
tx_genes <- tx_genes[o]

files <- list.files(pattern="*collapsed.bedgraph", recursive=TRUE)
start.time <- Sys.time()
for(i in files){
  dat <- read.table(i)
  #get identifier to name later file
  name= paste(head(unlist(strsplit(unlist(strsplit(i, "/"))[2], "_")), n=-1), collapse='_')
  print(name)
  #rename columns
  dat<- rename(dat,c("V1"="chromosome","V2"="start","V3"="count","V4"="x","V5"="strand"))
  #normalise counts based on total read count
  dat$count = dat$count/(sum(dat$count))*10000000
  dat$end <- as.integer(dat$start +1)
  #find Overlaps with genes by comparing genomic ranges
  dat_ranges <- makeGRangesFromDataFrame(dat, keep.extra.columns = TRUE)
  dat$overlaps=findOverlaps(dat_ranges, tx_genes, select = 'first')
  #lookup entrezId with overlap number
  dat$regions= with(dat, names(tx_genes)[overlaps])
  #set intergenic_$start/1000 instead of NA, to prevent loosing information
  dat[c("regions")][is.na(dat[c("regions")])] <- paste("intergenic",(dat[is.na(dat[c("regions")]),]$chromosome),(dat[is.na(dat[c("regions")]),]$strand),as.integer(dat[is.na(dat[c("regions")]),]$start/1000), sep="_")
  
  ###export dataframe from dat
  dat_export <- data.frame("geneID"=dat$regions, "Position"=dat$start, "BarcodeCount"=dat$count, "ReadCount"=dat$count, "PosFromAnno"=0, "Class"=0)
  dat_export <- dat_export[complete.cases(dat_export),]
  dat_export <- dat_export[with(dat_export, order(geneID, Position)),]
  #export data frame as csv, do not use quotes, rownames; do use colnames
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/data/Summary_counts/mouse/",name, ".summary.counts", sep=""), row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
  print(paste("/home/stroemic/hiwi_16/data/Summary_counts/mouse/",name, ".summary.counts", sep=""))
  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


