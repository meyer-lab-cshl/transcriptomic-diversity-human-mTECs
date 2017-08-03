library(data.table)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
library(dplyr)


my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chrM", "chrX", "chrY")
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
o <- order(as.numeric(tx_genes$gene_id))
tx_genes <- tx_genes[o]
transcript <- transcripts(txdb)

genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9


preprocess_feature_table <- function(datx, daty) {
  dat <- merge(datx, daty, by=c("geneID","Position"), all.x = TRUE )
  dat <- subset(dat, select = -c(ReadCount.x, PosFromAnno.x, Class.x, ReadCount.y, PosFromAnno.y, Class.y))
  dat$response<- "blub"
  dat[is.na(dat$BarcodeCount.y),]$response<- 'no'
  dat[!is.na(dat$BarcodeCount.y),]$response <- 'yes'
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
  dat <- dat[c(1,6,7,2,8,5,3)]
  return(dat)
  
}

internal <- read.table("/home/stroemic/hiwi_16/analysis/lars_perl/mouse/all_internal.isoforms.csv", sep=",", header = TRUE)
fantoms <- read.table("/home/stroemic/hiwi_16/analysis/lars_perl/mouse/all_fantoms.isoforms.csv", sep=",", header = TRUE)

total_int <- preprocess_feature_table(internal, fantoms)
total_int <- dplyr::rename(total_int, start=Position, Fantom5=response)
total_fa <- preprocess_feature_table(fantoms, internal)
total_fa <- dplyr::rename(total_fa, start=Position, Internal=response)


###Expressionlevel as rank for genes
genes_int <- total_int[,c("geneID","start","BarcodeCount.x")]
genes_int <- setDT(genes_int)[, lapply(.SD, sum), by=.(geneID), .SDcols=c("BarcodeCount.x")]
setDF(genes_int)
genes_int$rank <- rank(genes_int$BarcodeCount.x, ties.method = "average")
total_int$expr_rank_gene <- with(total_int, genes_int$rank[geneID])

###Expressionlevel relativ to gene and to total for each peak
total_int <- total_int %>% 
  group_by(geneID) %>% 
  mutate(sumCount = sum(BarcodeCount.x)) %>% 
  mutate(relCount = BarcodeCount.x/sumCount) %>% 
  group_by() %>% 
  mutate(relAllCount = as.numeric(BarcodeCount.x / sum(BarcodeCount.x))) %>%
  dplyr::select(-sumCount)


###Relative Position in gene
#GRanges(seqnames = test$chrom, strand = test$strand, ranges = IRanges( start = start(tx_genes[test$geneID]), end = test$Position ))

###Absolut Position
int_ranges <- makeGRangesFromDataFrame(total_int, keep.extra.columns = TRUE)

####Define GRanges
#TSS +-10
TSS10_genes <- promoters(tx_genes, upstream = 10, downstream = 10)
total_int$TSS_10_genes <- overlapsAny(int_ranges, TSS10_genes)
#TSS transcripts +-10
TSS10_transcripts <- promoters(transcript, upstream = 10, downstream = 10)
total_int$TSS_10_trans <- overlapsAny(int_ranges, TSS10_transcripts)
#TSS +-100
TSS100_genes <- promoters(tx_genes, upstream = 100, downstream = 100)
total_int$TSS_100_genes <- overlapsAny(int_ranges, TSS100_genes)
#TSS transcripts +-10
TSS100_transcripts <- promoters(transcript, upstream = 100, downstream = 100)
total_int$TSS_100_trans <- overlapsAny(int_ranges, TSS100_transcripts)

#First Exon 
exon = exonsBy(txdb, by='gene')
first_exon = lapply(exon, function(x) if( unique(strand(x)) == '+') x[1] else x[length(x)])
first_exon = do.call(GRangesList, first_exon)
first_exon = unlist(first_exon)
total_int$first_exon <- overlapsAny(int_ranges, first_exon)
#Other Exons 
other_exon = lapply(exon, function(x) if( unique(strand(x)) == '+') x[-1] else x[-length(x)])
other_exon = do.call(GRangesList, other_exon)
other_exon = unlist(other_exon)
total_int$other_exon <- overlapsAny(int_ranges, other_exon)
#Introns
all_exons <- unlist(exon)
introns <- GenomicRanges::setdiff(tx_genes, all_exons)
total_int$intron <- overlapsAny(int_ranges, introns)

#TTS +-100
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
plus_tts <- GRanges(seqnames = seqnames(plus), ranges=IRanges(start=(end(plus)-100), end =(end(plus)+100) ), strand = strand(plus))
minus_tts <- GRanges(seqnames = seqnames(minus), ranges=IRanges(start=(start(minus)-100), end =(start(minus)+100) ), strand = strand(minus))
TTS100_genes <- c(plus_tts, minus_tts)
total_int$TTS_100_genes <- overlapsAny(int_ranges, TTS100_genes)
#TTS +200
start(plus_tts) <- start(plus_tts) + 100
end(plus_tts) <- end(plus_tts) + 100
start(minus_tts) <- start(minus_tts) - 100
end(minus_tts) <- end(minus_tts) - 100
TTS200_genes <- c(plus_tts, minus_tts)
total_int$TTS_200_genes <- overlapsAny(int_ranges, TTS200_genes)

#TTS transcripts +-100
plus_trans <- transcript[strand(transcript) == '+']
minus_trans <- transcript[strand(transcript) == '-']
plus_trans_tts <- GRanges(seqnames = seqnames(plus_trans), ranges=IRanges(start=(end(plus_trans)-100), end =(end(plus_trans)+100) ), strand = strand(plus_trans))
minus_trans_tts <- GRanges(seqnames = seqnames(minus_trans), ranges=IRanges(start=(start(minus_trans)-100), end =(start(minus_trans)+100) ), strand = strand(minus_trans))
TTS100_trans <- c(plus_trans_tts, minus_trans_tts)
total_int$TTS_100_trans <- overlapsAny(int_ranges, TTS100_trans)
#TTS transcripts +200
start(plus_trans_tts) <- start(plus_trans_tts) + 100
end(plus_trans_tts) <- end(plus_trans_tts) + 100
start(minus_trans_tts) <- start(minus_trans_tts) - 100
end(minus_trans_tts) <- end(minus_trans_tts) - 100
TTS200_trans <- c(plus_trans_tts, minus_trans_tts)
total_int$TTS_200_trans <- overlapsAny(int_ranges, TTS200_trans)


#Downstream
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
start(plus) <- end(plus) 
end(plus) <- end(plus) + 1000
end(minus) <- start(minus)
start(minus) <- start(minus) - 1000
downstream_genes <- c(plus, minus)
total_int$downstream_gene <- overlapsAny(int_ranges, downstream_genes)
#Upstream
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
end(plus) <- start(plus) 
start(plus) <- start(plus) - 1000
start(minus) <- end(minus)
end(minus) <- end(minus) + 1000
upstream_genes <- c(plus, minus)
total_int$upstream_gene <- overlapsAny(int_ranges, upstream_genes)

#Antisense
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
strand(plus) <- '-'
strand(minus) <- '+'
antisense <- c(plus, minus)
total_int$antisense <- overlapsAny(int_ranges, antisense)

####2 Bases  after
total_int_plus <- total_int[total_int$strand == '+',]
total_int_minus <- total_int[total_int$strand == '-',]
total_int_plus$Basesafter <- getSeq(genome, total_int_plus$chrom, (total_int_plus$start+1), (total_int_plus$start +2), as.character=TRUE)
total_int_minus$Basesafter <- reverse(getSeq(genome, total_int_minus$chrom, (total_int_minus$start-2), (total_int_minus$start-1), as.character=TRUE))
total_int <- rbind(total_int_plus,total_int_minus)
####Base content 10bp window
total_int$Freq10A <- letterFrequency(getSeq(genome, total_int$chrom, (total_int$start-5), (total_int$start+5)), "A")
total_int$Freq10C <- letterFrequency(getSeq(genome, total_int$chrom, (total_int$start-5), (total_int$start+5)), "C")
total_int$Freq10G <- letterFrequency(getSeq(genome, total_int$chrom, (total_int$start-5), (total_int$start+5)), "G")
total_int$Freq10T <- letterFrequency(getSeq(genome, total_int$chrom, (total_int$start-5), (total_int$start+5)), "T")
####Base content 50bp window
total_int$Freq50A <- letterFrequency(getSeq(genome, total_int$chrom, (total_int$start-25), (total_int$start+25)), "A")
total_int$Freq50C <- letterFrequency(getSeq(genome, total_int$chrom, (total_int$start-25), (total_int$start+25)), "C")
total_int$Freq50G <- letterFrequency(getSeq(genome, total_int$chrom, (total_int$start-25), (total_int$start+25)), "G")
total_int$Freq50T <- letterFrequency(getSeq(genome, total_int$chrom, (total_int$start-25), (total_int$start+25)), "T")


#####Fantom5 table

###Expressionlevel as rank for genes
genes_fa <- total_fa[,c("geneID","start","BarcodeCount.x")]
genes_fa <- setDT(genes_fa)[, lapply(.SD, sum), by=.(geneID), .SDcols=c("BarcodeCount.x")]
setDF(genes_fa)
genes_fa$rank <- rank(genes_fa$BarcodeCount.x, ties.method = "average")
total_fa$expr_rank_gene <- with(total_fa, genes_fa$rank[geneID])

###Expressionlevel relativ to gene and to total for each peak
total_fa <- total_fa %>% 
  group_by(geneID) %>% 
  mutate(sumCount = sum(BarcodeCount.x)) %>% 
  mutate(relCount = BarcodeCount.x/sumCount) %>% 
  group_by() %>% 
  mutate(relAllCount = as.numeric(BarcodeCount.x / sum(BarcodeCount.x))) %>%
  dplyr::select(-sumCount)

###Absolut Position
fa_ranges <- makeGRangesFromDataFrame(total_fa, keep.extra.columns = TRUE)

####Define GRanges
#TSS +-10
total_fa$TSS_10_genes <- overlapsAny(fa_ranges, TSS10_genes)
#TSS transcripts +-10
total_fa$TSS_10_trans <- overlapsAny(fa_ranges, TSS10_transcripts)
#TSS +-100
total_fa$TSS_100_genes <- overlapsAny(fa_ranges, TSS100_genes)
#TSS transcripts +-10
total_fa$TSS_100_trans <- overlapsAny(fa_ranges, TSS100_transcripts)

#First Exon 
total_fa$first_exon <- overlapsAny(fa_ranges, first_exon)
#Other Exons 
total_fa$other_exon <- overlapsAny(fa_ranges, other_exon)
#Introns
total_fa$intron <- overlapsAny(fa_ranges, introns)

#TTS +-100
total_fa$TTS_100_genes <- overlapsAny(fa_ranges, TTS100_genes)
#TTS +200
total_fa$TTS_200_genes <- overlapsAny(fa_ranges, TTS200_genes)

#TTS transcripts +-100
total_fa$TTS_100_trans <- overlapsAny(fa_ranges, TTS100_trans)
#TTS transcripts +200
total_fa$TTS_200_trans <- overlapsAny(fa_ranges, TTS200_trans)


#Downstream
total_fa$downstream_gene <- overlapsAny(fa_ranges, downstream_genes)
#Upstream
total_fa$upstream_gene <- overlapsAny(fa_ranges, upstream_genes)

#Antisense
total_fa$antisense <- overlapsAny(fa_ranges, antisense)

####2 Bases  after
total_fa_plus <- total_fa[total_fa$strand == '+',]
total_fa_minus <- total_fa[total_fa$strand == '-',]
total_fa_plus$Basesafter <- getSeq(genome, total_fa_plus$chrom, (total_fa_plus$start+1), (total_fa_plus$start +2), as.character=TRUE)
total_fa_minus$Basesafter <- reverse(getSeq(genome, total_fa_minus$chrom, (total_fa_minus$start-2), (total_fa_minus$start-1), as.character=TRUE))
total_fa <- rbind(total_fa_plus,total_fa_minus)
####Base content 10bp window
total_fa$Freq10A <- letterFrequency(getSeq(genome, total_fa$chrom, (total_fa$start-5), (total_fa$start+5)), "A")
total_fa$Freq10C <- letterFrequency(getSeq(genome, total_fa$chrom, (total_fa$start-5), (total_fa$start+5)), "C")
total_fa$Freq10G <- letterFrequency(getSeq(genome, total_fa$chrom, (total_fa$start-5), (total_fa$start+5)), "G")
total_fa$Freq10T <- letterFrequency(getSeq(genome, total_fa$chrom, (total_fa$start-5), (total_fa$start+5)), "T")
####Base content 50bp window
total_fa <- total_fa[!(total_fa$start <25) ,]
total_fa$Freq50A <- letterFrequency(getSeq(genome, total_fa$chrom, (total_fa$start-25), (total_fa$start+25)), "A")
total_fa$Freq50C <- letterFrequency(getSeq(genome, total_fa$chrom, (total_fa$start-25), (total_fa$start+25)), "C")
total_fa$Freq50G <- letterFrequency(getSeq(genome, total_fa$chrom, (total_fa$start-25), (total_fa$start+25)), "G")
total_fa$Freq50T <- letterFrequency(getSeq(genome, total_fa$chrom, (total_fa$start-25), (total_fa$start+25)), "T")



write.table(total_int, file ="/home/stroemic/hiwi_16/analysis/r_random-forest/features_internal_data_mm9.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(total_fa, file ="/home/stroemic/hiwi_16/analysis/r_random-forest/features_fantom5_data_mm9.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")