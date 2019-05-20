library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)


setwd("/home/stroemic/hiwi_16/analysis/lars_perl/all_63")

#####Load Data structures
###annotation
my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
transcript <- transcripts(txdb)
###genome
genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
chromosome_length <- read.table("/home/stroemic/genomes/human/hg19/human.hg19.genome_length", stringsAsFactors=FALSE)
chromosome_length <- dplyr::rename(chromosome_length, chrom=V1, length=V2)
###files
files <- list.files(pattern = "isoforms.csv$")

###Define GRanges for Overlaps
#TSS +-10
TSS10_genes <- promoters(tx_genes, upstream = 10, downstream = 10)
#TSS +-100
TSS100_genes <- promoters(tx_genes, upstream = 100, downstream = 100)
#TSS transcripts +-10
TSS10_transcripts <- promoters(transcript, upstream = 10, downstream = 10)
#TSS transcripts +-10
TSS100_transcripts <- promoters(transcript, upstream = 100, downstream = 100)

#First Exon
exon = exonsBy(txdb, by='gene')
first_exon = lapply(exon, function(x) if( unique(strand(x)) == '+') x[1] else x[length(x)])
first_exon = do.call(GRangesList, first_exon)
first_exon = unlist(first_exon)
#Other Exons
other_exon = lapply(exon, function(x) if( unique(strand(x)) == '+') x[-1] else x[-length(x)])
other_exon = do.call(GRangesList, other_exon)
other_exon = unlist(other_exon)
#Introns
all_exons <- unlist(exon)
introns <- GenomicRanges::setdiff(tx_genes, all_exons)

#TTS +-100
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
plus_tts <- GRanges(seqnames = seqnames(plus), ranges=IRanges(start=(end(plus)-100), end =(end(plus)+100) ), strand = strand(plus))
minus_tts <- GRanges(seqnames = seqnames(minus), ranges=IRanges(start=(start(minus)-100), end =(start(minus)+100) ), strand = strand(minus))
TTS100_genes <- c(plus_tts, minus_tts)
#TTS +200
start(plus_tts) <- start(plus_tts) + 100
end(plus_tts) <- end(plus_tts) + 100
start(minus_tts) <- start(minus_tts) - 100
end(minus_tts) <- end(minus_tts) - 100
TTS200_genes <- c(plus_tts, minus_tts)
#TTS transcripts +-100
plus_trans <- transcript[strand(transcript) == '+']
minus_trans <- transcript[strand(transcript) == '-']
plus_trans_tts <- GRanges(seqnames = seqnames(plus_trans), ranges=IRanges(start=(end(plus_trans)-100), end =(end(plus_trans)+100) ), strand = strand(plus_trans))
minus_trans_tts <- GRanges(seqnames = seqnames(minus_trans), ranges=IRanges(start=(start(minus_trans)-100), end =(start(minus_trans)+100) ), strand = strand(minus_trans))
TTS100_trans <- c(plus_trans_tts, minus_trans_tts)
#TTS transcripts +200
start(plus_trans_tts) <- start(plus_trans_tts) + 100
end(plus_trans_tts) <- end(plus_trans_tts) + 100
start(minus_trans_tts) <- start(minus_trans_tts) - 100
end(minus_trans_tts) <- end(minus_trans_tts) - 100
TTS200_trans <- c(plus_trans_tts, minus_trans_tts)

#Downstream
start(plus) <- end(plus) 
end(plus) <- end(plus) + 1000
end(minus) <- start(minus)
start(minus) <- start(minus) - 1000
downstream_genes <- c(plus, minus)
#Upstream
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
end(plus) <- start(plus) 
start(plus) <- start(plus) - 1000
start(minus) <- end(minus)
end(minus) <- end(minus) + 1000
upstream_genes <- c(plus, minus)
#Antisense
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
strand(plus) <- '-'
strand(minus) <- '+'
antisense <- c(plus, minus)



preprocess_feature_table <- function(dat) {
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
  ###define end, adjust count based on strand
  dat$end <- as.integer(dat$Position +1)
  #dat[dat$strand == "-",]$count <- (dat[dat$strand == "-",]$count)*-1
  #rename columns and return new order
  dat <- dplyr::rename(dat, start=Position)
  dat <- dat[c(1,7,8,2,9,3)]
  return(dat)
  
}

start.time <- Sys.time()
for(i in files){
  dat <- read.table(i, header=TRUE, sep = ",")
  #get identifier to name later file
  if (grepl("_", i)) {
    name= paste((head(unlist(strsplit(i, "[_]" )), n=-1)), collapse='_')
    print(name)
  } else {
    name= paste((head(unlist(strsplit(i, "[.]" )), n=-2)), collapse='.')
    print(name)
  }
  dat <- preprocess_feature_table(dat)
  ##Define features
  #Expressionlevel as rank for genes
  genes_dat <- dat[,c("geneID","start","BarcodeCount")]
  genes_dat <- setDT(genes_dat)[, lapply(.SD, sum), by=.(geneID), .SDcols=c("BarcodeCount")]
  setDF(genes_dat)
  genes_dat$rank <- rank(genes_dat$BarcodeCount, ties.method = "average")
  dat$expr_rank_gene <- with(dat, genes_dat$rank[geneID])
  #Expressionlevel relativ to gene and to total for each peak
  dat <- dat %>%
    group_by(geneID) %>%
    mutate(sumCount = sum(BarcodeCount)) %>%
    mutate(relCount = BarcodeCount/sumCount) %>%
    group_by() %>%
    mutate(relAllCount = as.numeric(BarcodeCount / sum(BarcodeCount))) %>%
    dplyr::select(-sumCount)
  ###Absolut Positios
  dat_ranges <- makeGRangesFromDataFrame(dat, keep.extra.columns = TRUE)
  dat$TSS_10_genes <- overlapsAny(dat_ranges, TSS10_genes)
  dat$TSS_100_genes <- overlapsAny(dat_ranges, TSS100_genes)
  dat$TSS_10_trans <- overlapsAny(dat_ranges, TSS10_transcripts)
  dat$TSS_100_trans <- overlapsAny(dat_ranges, TSS100_transcripts)
  dat$first_exon <- overlapsAny(dat_ranges, first_exon)
  dat$other_exon <- overlapsAny(dat_ranges, other_exon)
  dat$intron <- overlapsAny(dat_ranges, introns)
  dat$TTS_100_genes <- overlapsAny(dat_ranges, TTS100_genes)
  dat$TTS_200_genes <- overlapsAny(dat_ranges, TTS200_genes)
  dat$TTS_100_trans <- overlapsAny(dat_ranges, TTS100_trans)
  dat$TTS_200_trans <- overlapsAny(dat_ranges, TTS200_trans)
  dat$downstream_gene <- overlapsAny(dat_ranges, downstream_genes)
  dat$upstream_gene <- overlapsAny(dat_ranges, upstream_genes)
  dat$antisense <- overlapsAny(dat_ranges, antisense)
  
  ####3 Bases and after
  print(dim(dat))
  dat <- merge(dat, chromosome_length,by="chrom" ,all.x =TRUE)
  dat <- dat[(dat$start > 25) & (dat$start < (dat$length - 24)),]
  print(dim(dat))
  dat_plus <- dat[dat$strand == '+',]
  dat_minus <- dat[dat$strand == '-',]
  dat_plus$Basesafter <- getSeq(genome, dat_plus$chrom, (dat_plus$start+1), (dat_plus$start +2), as.character=TRUE)
  dat_minus$Basesafter <- reverse(getSeq(genome, dat_minus$chrom, (dat_minus$start-2), (dat_minus$start-1), as.character=TRUE))
  dat <- rbind(dat_plus,dat_minus)
  ####Base content 10bp window
  dat$wind10 <- getSeq(genome, dat$chrom, (dat$start-5), (dat$start+5))
  dat$Freq10A <- letterFrequency(dat$wind10, "A")
  dat$Freq10C <- letterFrequency(dat$wind10, "C")
  dat$Freq10G <- letterFrequency(dat$wind10, "G")
  dat$Freq10T <- letterFrequency(dat$wind10, "T")
  ####Base content 50bp window
  dat$wind50 <- getSeq(genome, dat$chrom, (dat$start-25), (dat$start+25))
  dat$Freq50A <- letterFrequency(dat$wind50, "A")
  dat$Freq50C <- letterFrequency(dat$wind50, "C")
  dat$Freq50G <- letterFrequency(dat$wind50, "G")
  dat$Freq50T <- letterFrequency(dat$wind50, "T")
  dat <- subset(dat, select= -c(wind10, wind50, length))
  
  #export data frame as csv, do not use quotes, rownames; do use colnames
  write.table(dat, file =paste("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/features/",name, ".features.csv", sep=""), row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
  print(paste("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/features/",name, ".features.csv", sep=""))
  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
