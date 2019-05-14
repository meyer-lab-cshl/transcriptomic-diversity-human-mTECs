library(dplyr)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

#### functions ####
grange_to_df <- function(x) {
  name <- deparse(substitute(x))
  df <- data.frame(seqnames=seqnames(x),
             starts=start(x)-1,
             ends=end(x),
             names=c(rep(".", length(x))),
             scores=c(rep(".", length(x))),
             strands=strand(x))
  df <- df[with(df, order(seqnames, starts)),]
  write.table(df, file=paste("/home/stroemic/hiwi_16/data/bed_for_igv/",name,".bed", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
  return(df)
}

#### read in data ####
##### BioMart object
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="grch37.ensembl.org", 
                   dataset="hsapiens_gene_ensembl")

##### txdb gene object
my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
### keep only protein coding genes 
mart_gene_type <- getBM(attributes = c('entrezgene', 'gene_biotype'), 
                        filters = 'entrezgene',
                        values = names(tx_genes), 
                        mart=ensembl)
mart_gene_type <- mart_gene_type[grepl("protein_coding", mart_gene_type$gene_biotype, ignore.case = TRUE),]
tx_genes <- tx_genes[(elementMetadata(tx_genes)$gene_id %in% mart_gene_type$entrezgene)]
### sort by geneID
o <- order(as.numeric(tx_genes$gene_id))
tx_genes <- tx_genes[o]

#### Define GRanges ####
#TSS +-10
TSS10_genes <- promoters(tx_genes, upstream = 10, downstream = 10)
#TSS +-100
TSS100_genes <- promoters(tx_genes, upstream = 100, downstream = 100)


#First Exon
exon = exonsBy(txdb, by='gene')
exon <- exon[names(exon) %in% mart_gene_type$entrezgene]
oe <- order(as.numeric(names(exon)))
exon <- exon[oe]
first_exon = lapply(exon, function(x) if( unique(strand(x)) == '+') x[1] else x[length(x)])
first_exon = do.call(GRangesList, first_exon)
first_exon = unlist(first_exon)
#Other Exons
other_exon = lapply(exon, function(x) if( unique(strand(x)) == '+') x[-1] else x[-length(x)])
other_exon = do.call(GRangesList, other_exon)
other_exon = unlist(other_exon)

#TSS10 & 100 other exon
TSS10_other_exon <- promoters(other_exon, upstream = 10, downstream = 10)
TSS100_other_exon <- promoters(other_exon, upstream = 100, downstream = 100)

#Introns
# introns <- sapply(seq_along(exon), function(i) GenomicRanges::setdiff(tx_genes[names(exon[i])], exon[[i]]), USE.NAMES = TRUE)
# names(introns) <- names(exon)
# introns <- do.call(GRangesList, introns)
# introns = unlist(introns)
introns <- readRDS(file="/home/stroemic/hiwi_16/analysis/gene_lists/introns.rds")

#TTS +-100
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
start(plus) <- end(plus) -100
end(plus) <- end(plus) + 100
end(minus) <- start(minus) +100
start(minus) <- start(minus) - 100
TTS100_genes <- c(plus, minus)
otts <- order(as.numeric(names(TTS100_genes)))
TTS100_genes <- TTS100_genes[otts]

#Downstream
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
start(plus) <- end(plus) 
end(plus) <- end(plus) + 1000
end(minus) <- start(minus)
start(minus) <- start(minus) - 1000
downstream_genes <- c(plus, minus)
od <- order(as.numeric(names(downstream_genes)))
downstream_genes <- downstream_genes[od]
#Upstream
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
end(plus) <- start(plus) 
start(plus) <- start(plus) - 1000
start(minus) <- end(minus)
end(minus) <- end(minus) + 1000
upstream_genes <- c(plus, minus)
ou <- order(as.numeric(names(upstream_genes)))
upstream_genes <- upstream_genes[ou]
#Antisense
plus <- tx_genes[strand(tx_genes) == '+']
minus <- tx_genes[strand(tx_genes) == '-']
strand(plus) <- '-'
strand(minus) <- '+'
antisense <- c(plus, minus)
antisense <- antisense + 1000
oa <- order(as.numeric(names(antisense)))
antisense <- antisense[oa]


#### analysis ####
# write granges to dfs and export as bed files 
df_tss10 <- grange_to_df(TSS10_genes)
df_tss100 <- grange_to_df(TSS100_genes)
df_first_exon <- grange_to_df(first_exon)
other_exon_red <- reduce(other_exon)
df_other_exon_red <- grange_to_df(other_exon_red)
df_tss10_oe <- grange_to_df(TSS10_other_exon)
df_tss100_oe <- grange_to_df(TSS100_other_exon)
df_introns <- grange_to_df(introns)
df_tts100 <- grange_to_df(TTS100_genes)
df_downstream <- grange_to_df(downstream_genes)
df_upstream <- grange_to_df(upstream_genes)
df_antisense <- grange_to_df(antisense)