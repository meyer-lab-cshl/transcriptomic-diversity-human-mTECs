library(plyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

#####Load gene object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)

#####consensus sites csv read in
dat <- read.table("/home/stroemic/hiwi_16/analysis/lars_perl/all_63_consensus_sites.csv", sep = ",")
dat <- rename(dat,c("V1"="geneID","V2"="start","V3"="count"))
dat <- subset(dat, select = -c(V4, V5, V6))
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
#####export dataframe from dat
dat_export <- data.frame("chrom"=dat$chrom, "start"=dat$start, "end"=dat$end, "dataValue"=dat$count)
dat_export <- dat_export[complete.cases(dat_export),]
write.table(dat_export, file ="/home/stroemic/hiwi_16/analysis/lars_perl/all_63_consensus_sites.bedgraph", row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
