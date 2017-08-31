library(plyr)
library(data.table)
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
datalist <- list()
c = 1
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
  dat_agg <- setDT(dat)[, lapply(.SD, sum), by=.(chr,strand,pos), .SDcols=c("count")]
  setDF(dat_agg)
  dat_agg <- plyr::rename(dat_agg, c("pos"="start"))
  #dat<- aggregate(dat$count, by=list(dat$chr, dat$pos, dat$strand), FUN=sum)
  #rename columns
  dat_agg$end <- as.integer(dat_agg$start +1)
  #normalise counts based on total read count
  dat_agg$count = dat_agg$count/(sum(dat_agg$count))*10000000
  #find Overlaps with genes by comparing genomic ranges
  dat_ranges <- makeGRangesFromDataFrame(dat_agg[,c(1,2,3,5,4)], keep.extra.columns = TRUE)
  dat_agg$overlaps=findOverlaps(dat_ranges, tx_genes, select = 'first')
  #lookup entrezId with overlap number
  dat_agg$regions= with(dat_agg, names(tx_genes)[overlaps])
  #set intergenic_$start/1000 instead of NA, to prevent loosing information
  dat_agg[c("regions")][is.na(dat_agg[c("regions")])] <- paste("intergenic",(dat_agg[is.na(dat_agg[c("regions")]),]$chr),(dat_agg[is.na(dat_agg[c("regions")]),]$strand),as.integer(dat_agg[is.na(dat_agg[c("regions")]),]$start/1000), sep="_")
  ###export dataframe from dat
  dat_export <- data.frame("chrom"=dat_agg$chr, "geneID"=dat_agg$regions, "strand"=dat_agg$strand, "start"=dat_agg$start, "end"=dat_agg$end, "BarcodeCount"=dat_agg$count)
  dat_export <- dat_export[complete.cases(dat_export),]
  dat_export <- dat_export[with(dat_export, order(geneID, start)),]
  #export data frame as csv, do not use quotes, rownames; do use colnames
  write.table(dat_export, file =paste("/home/stroemic/hiwi_16/data/raw_positions/",name, ".positions.csv", sep=""), row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

  dat_export$iteration <- c
  datalist[[c]] <- dat_export
  c <- c+1

}


 print("rbinding")
 dat_all = do.call(rbind, datalist)
 
 print("aggregating")
 dat_agg_all <- setDT(dat_all)[, lapply(.SD, sum), by=.(chrom,geneID,strand,start,end), .SDcols=c("BarcodeCount")]
 setDF(dat_agg_all)
 
 print("writing out")
 write.table(dat_agg_all, file ="/home/stroemic/hiwi_16/data/raw_positions/all_mTECs.positions.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
