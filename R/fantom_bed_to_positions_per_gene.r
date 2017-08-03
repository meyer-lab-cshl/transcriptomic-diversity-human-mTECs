library(plyr)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

setwd("/home/stroemic/hiwi_16/data/external/Fantom5/CTSS_BED/tissues/")

my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

#####Load gene object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
tx_genes <- tx_genes + 100
o <- order(as.numeric(tx_genes$gene_id))
tx_genes <- tx_genes[o]

fantoms <- list.files(pattern="*.bed.gz")

datalist<-list()
datalist_wo<-list()
c <- 1

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
  
  ###export dataframe from dat chrom,geneID,strand,start,end,BarcodeCount

  dat_export <- data.frame("chrom"=dat$chromosome, "geneID"=dat$regions, "strand"=dat$strand, "start"=dat$start, "end"=dat$end, "BarcodeCount"=dat$count)
  dat_export <- dat_export[complete.cases(dat_export),]
  dat_export <- dat_export[with(dat_export, order(chrom, geneID)),]
  #export data frame as csv, do not use quotes, rownames; do use colnames
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/data/raw_positions/",name, ".positions.csv", sep=""), row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
  
  dat_export$iteration <- c
  if (grepl("thymus", i)) {
    print("one")
    datalist[[c]] <- dat_export
  } else {
    print("both")
    datalist[[c]] <- dat_export
    datalist_wo[[c]] <- dat_export
  }
  c <- c+1
}

print("rbinding")
dat_fantom5 = do.call(rbind, datalist)
dat_fantom5_wo_thymus = do.call(rbind, datalist_wo)

print("aggregating")
dat_agg_fantom5 <- setDT(dat_fantom5)[, lapply(.SD, sum), by=.(chrom,geneID,strand,start,end), .SDcols=c("BarcodeCount")]
setDF(dat_agg_fantom5)
dat_agg_fantom5_wo_thymus <- setDT(dat_fantom5_wo_thymus)[, lapply(.SD, sum), by=.(chrom,geneID,strand,start,end), .SDcols=c("BarcodeCount")]
setDF(dat_agg_fantom5_wo_thymus)

print("writing out")
write.table(dat_agg_fantom5, file ="/home/stroemic/hiwi_16/data/raw_positions/all.tissues.positions.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(dat_agg_fantom5_wo_thymus, file ="/home/stroemic/hiwi_16/data/raw_positions/all.tissues.wo.thymus.positions.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

