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

filter_df <- function(x, window_size) {
  filt <- x %>% group_by(geneID) %>% mutate(the_rank  = rank(-BarcodeCount, ties.method = "first")) %>% filter(the_rank == 1)
  filt_plus <- filt[filt$strand == '+',]
  filt_minus <- filt[filt$strand == '-',]
  filt_plus$sequence <- getSeq(genome, filt_plus$chrom, (filt_plus$start-window_size), (filt_plus$start +window_size), as.character=TRUE)
  filt_plus$ann_sequence <- getSeq(genome, as.character(seqnames(tx_genes[filt_plus$geneID])), (start(tx_genes[filt_plus$geneID])-window_size), (start(tx_genes[filt_plus$geneID]) +window_size), as.character=TRUE )
  filt_minus$sequence <- reverse(getSeq(genome, filt_minus$chrom, (filt_minus$start-window_size), (filt_minus$start+window_size), as.character=TRUE))
  filt_minus$ann_sequence <- reverse(getSeq(genome, as.character(seqnames(tx_genes[filt_minus$geneID])), (end(tx_genes[filt_minus$geneID])-window_size), (end(tx_genes[filt_minus$geneID]) +window_size), as.character=TRUE ))
  
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
##table with housekeeping genes defined on all tissues
hk_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/housekeeping_genes.csv", header=TRUE, sep=",")

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
mart_gene_type <- mart_gene_type[grepl("protein_coding", mart_gene_type$gene_biotype, ignore.case = TRUE),]
tx_genes <- tx_genes[(elementMetadata(tx_genes)$gene_id %in% mart_gene_type$entrezgene)]
### sort by geneID
o <- order(as.numeric(tx_genes$gene_id))
tx_genes <- tx_genes[o]


#### Define gene classes ####
aire_dep_tras <- merge(tra_genes, aire_dep_genes, by="entrezgene")
fezf2_dep_tras <- merge(tra_genes, fezf2_dep_genes, by="entrezgene")

# other tras
other_tras <- merge(tra_genes, aire_dep_genes, by="entrezgene", all.x = TRUE)
other_tras <- other_tras[is.na(other_tras$ensembl_gene_id.y),]
other_tras <- merge(other_tras, fezf2_dep_genes, ny="entrezgene", all.x = TRUE)
other_tras <- other_tras[is.na(other_tras$ensembl_gene_id),]

#other genes
others <- names(tx_genes)[(!names(tx_genes) %in% tra_genes$entrezgene)]
others <- others[!(others %in% hk_genes$entrezgene)]
others <- as.integer(others)


##### human 
###filtered
##mtec all positions
mTEC_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered.csv", header=TRUE, sep=",")
##tissue (wo thymus) all positions
tissues_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.wo.thymus.filtered.csv", header = TRUE, sep = ",")



#### analysis ####
genes_entrez_mTECs <- list(aire_dep_tras$entrezgene, fezf2_dep_tras$entrezgene, other_tras$entrezgene, hk_genes$entrezgene , others)
genes_entrez_tissues <- list(aire_dep_tras$entrezgene, fezf2_dep_tras$entrezgene, other_tras$entrezgene, hk_genes$entrezgene, others)
 
gene_classes_mTECs <- c("aire_dep_tras_mTECs", "fezf2_dep_tras_mTECs", "other_tras_mTECs", "hk_mTECs", "other_genes_mTECs")
gene_classes_tissues <- c("aire_dep_tras_tissues", "fezf2_dep_tras_tissues", "other_tras_tissues", "hk_tissues"  ,"other_genes_tissues")


mTEC_pos_fil <- preprocess_df(mTEC_pos_fil)
tissues_pos_fil <- preprocess_df(tissues_pos_fil)

filtered_pos_mTECs_51 <- filter_df(mTEC_pos_fil, 25)
filtered_pos_tissues_51 <- filter_df(tissues_pos_fil, 25)

filtered_pos_mTECs_11 <- filter_df(mTEC_pos_fil, 5)
filtered_pos_tissues_11 <- filter_df(tissues_pos_fil, 5)


##window of 51
for (i in 1:5) {
  name <- gene_classes_mTECs[[i]]
  print(name)
  df_main <- as.data.frame(filtered_pos_mTECs_51[filtered_pos_mTECs_51$geneID %in% genes_entrez_mTECs[[i]],c("name","sequence")])
  df_anno <- as.data.frame(filtered_pos_mTECs_51[filtered_pos_mTECs_51$geneID %in% genes_entrez_mTECs[[i]],c("name","ann_sequence")])
  print(dim(df_main))
  dd <- dataframe2fas(df_main, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences_mainTSS_51.fa", sep=""))
  dd <- dataframe2fas(df_anno, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences_annotatedTSS_51.fa", sep=""))
  
}

for (i in 1:5) {
  name <- gene_classes_tissues[[i]]
  print(name)
  df_main <- as.data.frame(filtered_pos_tissues_51[filtered_pos_tissues_51$geneID %in% genes_entrez_tissues[[i]],c("name","sequence")])
  df_anno <- as.data.frame(filtered_pos_tissues_51[filtered_pos_tissues_51$geneID %in% genes_entrez_tissues[[i]],c("name","ann_sequence")])
  print(dim(df_main))
  dd <- dataframe2fas(df_main, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences_mainTSS_51.fa", sep=""))
  dd <- dataframe2fas(df_anno, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences_annotatedTSS_51.fa", sep=""))
}

## window of 11
for (i in 1:5) {
  name <- gene_classes_mTECs[[i]]
  print(name)
  df_main <- as.data.frame(filtered_pos_mTECs_11[filtered_pos_mTECs_11$geneID %in% genes_entrez_mTECs[[i]],c("name","sequence")])
  df_anno <- as.data.frame(filtered_pos_mTECs_11[filtered_pos_mTECs_11$geneID %in% genes_entrez_mTECs[[i]],c("name","ann_sequence")])
  print(dim(df_main))
  dd <- dataframe2fas(df_main, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences_mainTSS_11.fa", sep=""))
  dd <- dataframe2fas(df_anno, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences_annotatedTSS_11.fa", sep=""))
  
}

for (i in 1:5) {
  name <- gene_classes_tissues[[i]]
  print(name)
  df_main <- as.data.frame(filtered_pos_tissues_11[filtered_pos_tissues_11$geneID %in% genes_entrez_tissues[[i]],c("name","sequence")])
  df_anno <- as.data.frame(filtered_pos_tissues_11[filtered_pos_tissues_11$geneID %in% genes_entrez_tissues[[i]],c("name","ann_sequence")])
  print(dim(df_main))
  dd <- dataframe2fas(df_main, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences_mainTSS_11.fa", sep=""))
  dd <- dataframe2fas(df_anno, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/sequences/",name,"_sequences_annotatedTSS_11.fa", sep=""))
}

df_mTECs_51_main <- as.data.frame(filtered_pos_mTECs_51[c("name","sequence")])
df_mTECs_51_anno <- as.data.frame(filtered_pos_mTECs_51[c("name","ann_sequence")])
df_mTECs_11_main <- as.data.frame(filtered_pos_mTECs_11[c("name","sequence")])
df_mTECs_11_anno <- as.data.frame(filtered_pos_mTECs_11[c("name","ann_sequence")])

temp1 <- dataframe2fas(df_mTECs_51_main, file="/home/stroemic/hiwi_16/analysis/gene_lists/sequences/all_mTEC_protein_coding_genes_sequences_mainTSS_51.fa")
temp2 <- dataframe2fas(df_mTECs_51_anno, file="/home/stroemic/hiwi_16/analysis/gene_lists/sequences/all_mTEC_protein_coding_genes_sequences_annotatedTSS_51.fa")
temp3 <- dataframe2fas(df_mTECs_11_main, file="/home/stroemic/hiwi_16/analysis/gene_lists/sequences/all_mTEC_protein_coding_genes_sequences_mainTSS_11.fa")
temp4 <- dataframe2fas(df_mTECs_11_anno, file="/home/stroemic/hiwi_16/analysis/gene_lists/sequences/all_mTEC_protein_coding_genes_sequences_annotatedTSS_11.fa")

