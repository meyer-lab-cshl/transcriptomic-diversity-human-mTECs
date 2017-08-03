library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqRFLP)
library(biomaRt)

#### functions ####
preprocess_df <- function(x) {
  x <- x[x$geneID %in% names(tx_genes),]
  x <- droplevels(x)
  x$geneID <- as.character(x$geneID)
  
  x_plus <- x[x$strand == "+",]
  x_minus <- x[x$strand == "-",]
  
  x_plus$rel_pos <- with(x_plus, (x_plus$start - start(tx_genes[x_plus$geneID])))
  x_minus$rel_pos <- with(x_minus, ( end(tx_genes[x_minus$geneID]) - x_minus$start))
  
  x <- rbind(x_plus, x_minus)
  x <- x[with(x, order(chrom, geneID, rel_pos)),]
  x$rel_pos <- as.character(x$rel_pos)
  x$name <- paste(x$chrom, x$geneID, x$strand, x$start, sep = "_")
  return(x)
}

filter_df <- function(x) {
  filt <- x %>% group_by(geneID) %>% mutate(the_rank  = rank(-BarcodeCount, ties.method = "first")) %>% filter(the_rank == 1)
  
  filt_plus <- filt[filt$strand == '+',]
  filt_minus <- filt[filt$strand == '-',]
  filt_plus$sequence <- getSeq(genome, filt_plus$chrom, (filt_plus$start-50), (filt_plus$start +50), as.character=TRUE)
  filt_minus$sequence <- reverse(getSeq(genome, filt_minus$chrom, (filt_minus$start-50), (filt_minus$start+50), as.character=TRUE))
  filt <- rbind(filt_plus,filt_minus)
  
  return(filt)
}


#### read in data ####
##table with tras
tra_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_tra_genes_sansom.csv", header=TRUE, sep=",")
##table with aire dep genes 
aire_dep_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_aire_dep_san_genes.csv", header=TRUE, sep=",")
##table with fezf2 dep genes 
fezf2_dep_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_fezf2_dep_genes.csv", header=TRUE, sep=",")
##table with thymus specific genes 
thymus_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/thymus_specific_genes.csv", header=TRUE, sep=",")
thymus_genes <- data.frame("entrezgene" =thymus_genes[!grepl("inter", thymus_genes$entrezgene),])
##table with mTEC specific genes
mTEC_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/mTEC_specific_genes.csv", header=TRUE, sep=",")
mTEC_genes <- data.frame("entrezgene" = mTEC_genes[!grepl("inter", mTEC_genes$entrezgene),])
##table with housekeeping genes defined on all tissues
hk_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/housekeeping_genes.csv", header=TRUE, sep=",")

### Define gene classes ###
aire_dep_tras <- merge(tra_genes, aire_dep_genes, by="entrezgene")
fezf2_dep_tras <- merge(tra_genes, fezf2_dep_genes, by="entrezgene")

# other tras
other_tras <- merge(tra_genes, aire_dep_genes, by="entrezgene", all.x = TRUE)
other_tras <- other_tras[is.na(other_tras$ensembl_gene_id.y),]
other_tras <- merge(other_tras, fezf2_dep_genes, ny="entrezgene", all.x = TRUE)
other_tras <- other_tras[is.na(other_tras$ensembl_gene_id),]


##### human 
###filtered
##mtec all positions
mTEC_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered.csv", header=TRUE, sep=",")
##tissue (wo thymus) all positions
tissues_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.wo.thymus.filtered.csv", header = TRUE, sep = ",")
###unfiltered
##mtec features 
mTEC_features_unfil <- read.table("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/features/all_mTECs.features.csv", header=TRUE, sep=",")
mTEC_pos_unfil <- mTEC_features_unfil[,c(1:6)]
##tissue (wo thymus) features
tissues_features_unfil <- read.table("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/features/all.tissues.wo.thymus.features.csv", header = TRUE, sep = ",")
tissues_pos_unfil <- tissues_features_unfil[,c(1:6)]


features_dfs <- list(human_filt=list(mTEC_pos_fil, tissues_pos_fil), 
                     human_unfilt=list(mTEC_pos_unfil, tissues_pos_unfil))


#### build genomes, txdb and biomart ####
##### genome object
genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

##### BioMart object
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="grch37.ensembl.org", 
                   dataset="hsapiens_gene_ensembl")

##### txdb gene object
my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
### filter out non-coding RNA (snoRNA, lincRNA, miRNA, snRNA, misc_RNA, 3prime_overlapping_ncrna)
mart_gene_type <- getBM(attributes = c('entrezgene', 'gene_biotype'), 
                        filters = 'entrezgene',
                        values = names(tx_genes), 
                        mart=ensembl)
mart_gene_type <- mart_gene_type[!grepl("rna", mart_gene_type$gene_biotype, ignore.case = TRUE),]
tx_genes <- tx_genes[(elementMetadata(tx_genes)$gene_id %in% mart_gene_type$entrezgene)]
### sort by geneID
o <- order(as.numeric(tx_genes$gene_id))
tx_genes <- tx_genes[o]


#### analysis ####
genes_entrez_mTECs <- list(aire_dep_tras$entrezgene, fezf2_dep_tras$entrezgene, other_tras$entrezgene, hk_genes$entrezgene ,mTEC_genes$entrezgene, thymus_genes$entrezgene)
genes_entrez_tissues <- list(hk_genes$entrezgene, aire_dep_tras$entrezgene, fezf2_dep_tras$entrezgene, other_tras$entrezgene)
 
gene_classes_mTECs <- c("aire_dep_tras_mTECs", "fezf2_dep_tras_mTECs", "other_tras_mTECs", "hk_mTECs", "mTEC_spec_mTECs", "thymus_spec_mTECs")
gene_classes_tissues <- c("hk_tissues", "aire_dep_tras_tissues", "fezf2_dep_tras_tissues", "other_tras_tissues")


mTEC_pos_fil <- preprocess_df(mTEC_pos_fil)
tissues_pos_fil <- preprocess_df(tissues_pos_fil)

filtered_pos_mTECs <- filter_df(mTEC_pos_fil)
filtered_pos_tissues <- filter_df(tissues_pos_fil)


for (i in 1:6) {
  name <- gene_classes_mTECs[[i]]
  print(name)
  df <- as.data.frame(filtered_pos_mTECs[filtered_pos_mTECs$geneID %in% genes_entrez_mTECs[[i]],c("name","sequence")])
  print(dim(df))
  dd <- dataframe2fas(df, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences.fa", sep=""))
}

for (i in 1:4) {
  name <- gene_classes_tissues[[i]]
  print(name)
  df <- as.data.frame(filtered_pos_tissues[filtered_pos_tissues$geneID %in% genes_entrez_tissues[[i]],c("name","sequence")])
  print(dim(df))
  dd <- dataframe2fas(df, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences.fa", sep=""))
}
