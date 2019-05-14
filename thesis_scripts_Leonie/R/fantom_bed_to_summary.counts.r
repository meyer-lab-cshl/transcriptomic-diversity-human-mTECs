library(plyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

setwd("/home/stroemic/hiwi_16/data/external/Fantom5/CTSS_BED/tissues/")

my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

#####Load gene object order by entrez gene ID to make sure findOverlaps= 'first' takes smaller ID
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
tx_genes <- tx_genes + 100
o <- order(as.numeric(tx_genes$gene_id))
tx_genes <- tx_genes[o]

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
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/data/Summary_counts/",name, ".summary.counts", sep=""), row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
  print(paste("/home/stroemic/hiwi_16/data/Summary_counts/",name, ".summary.counts", sep=""))
   
}
