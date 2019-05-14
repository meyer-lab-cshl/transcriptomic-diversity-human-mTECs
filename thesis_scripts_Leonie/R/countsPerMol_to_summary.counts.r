library(plyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

setwd("/home/stroemic/hiwi_16/analysis/STAR_alignment_hg19/5Seq_processed_IVT/")

my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

#####Load gene object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
tx_genes <- tx_genes + 100
o <- order(as.numeric(tx_genes$gene_id))
tx_genes <- tx_genes[o]

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