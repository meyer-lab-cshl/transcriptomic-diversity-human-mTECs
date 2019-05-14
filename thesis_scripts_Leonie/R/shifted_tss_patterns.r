library(dplyr)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)
library(reshape2)

#### functions ####
preprocess_df <- function(x, keep) {
  x <- x[x$geneID %in% names(cds),]
  x <- x[x$geneID %in% keep,]
  x <- droplevels(x)
  x$geneID <- as.character(x$geneID)
  return(x)
}

#### read in ####
#### data tables
tra_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_tra_genes_sansom.csv", header=TRUE, sep=",")
##mtec all positions
mTEC_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered.csv", header=TRUE, sep=",")
##tissue (wo thymus) all positions
tissues_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.wo.thymus.filtered.csv", header = TRUE, sep = ",")
##get common genes for calculation
keep <- intersect(mTEC_pos_fil$geneID, tissues_pos_fil$geneID)

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

##### txdb cds object
cds <- cdsBy(txdb, by="gene")
cds <- cds[names(cds) %in% mart_gene_type$entrezgene]
cds <- cds[names(cds) %in% keep]
oc <- order(as.numeric(names(cds)))
cds <- cds[oc]
first_cds = lapply(cds, function(x) if( unique(strand(x)) == '+') x[1] else x[length(x)])
first_cds = do.call(GRangesList, first_cds)
first_cds = unlist(first_cds)
last_cds = lapply(cds, function(x) if( unique(strand(x)) == '+') x[length(x)] else x[1])
last_cds = do.call(GRangesList, last_cds)
last_cds = unlist(last_cds)


#### analysis ####

mTEC_pos_fil <- preprocess_df(mTEC_pos_fil, keep)
tissues_pos_fil <- preprocess_df(tissues_pos_fil, keep)

list_genes_names <- names(first_cds)
list_genes <- vector("list", length(list_genes_names))
names(list_genes) <- list_genes_names

start.time <- Sys.time()
for (i in names(list_genes)) {
  print(i)
  #build matrix to write counts into 
  matrix_gene <- matrix(0, nrow=2, ncol=3)
  rownames(matrix_gene) <-c("tissues","mTEC") 
  colnames(matrix_gene) <-c("upstream","cds","downstream")
  #subdata frames for the gene
  mTEC_gene <- mTEC_pos_fil[mTEC_pos_fil$geneID == i,]
  tissues_gene <- tissues_pos_fil[tissues_pos_fil$geneID == i,]
  #iterate over all positions in mTEC data for the gene 
  for (pos_i in 1:nrow(mTEC_gene)) {
    #### keep in mind gene sense or antisense
    if (unique(strand(first_cds[i])) == '+'){
      #start and end of cds from cds objects
      start_cds <- start(first_cds[i])
      end_cds <- end(last_cds[i])
      #### counts in matrix according to position either in cds, upstream or downstream
      if ((mTEC_gene$start[pos_i] >= start_cds) & (mTEC_gene$start[pos_i] <= end_cds)) {
        matrix_gene["mTEC","cds"] <- matrix_gene["mTEC","cds"] + mTEC_gene$BarcodeCount[pos_i]
      } else if (mTEC_gene$start[pos_i] < start_cds) {
        matrix_gene["mTEC","upstream"] <- matrix_gene["mTEC","upstream"] + mTEC_gene$BarcodeCount[pos_i]
      } else {
        matrix_gene["mTEC","downstream"] <- matrix_gene["mTEC","downstream"] + mTEC_gene$BarcodeCount[pos_i]
      }
      
    } else {
      ## negative strand with different start and end 
      start_cds <- end(first_cds[i])
      end_cds <- start(last_cds[i])
      if ((mTEC_gene$start[pos_i] <= start_cds) & (mTEC_gene$start[pos_i] >= end_cds)) {
        matrix_gene["mTEC","cds"] <- matrix_gene["mTEC","cds"] + mTEC_gene$BarcodeCount[pos_i]
      } else if (mTEC_gene$start[pos_i] > start_cds) {
        matrix_gene["mTEC","upstream"] <- matrix_gene["mTEC","upstream"] + mTEC_gene$BarcodeCount[pos_i]
      } else { 
        matrix_gene["mTEC","downstream"] <- matrix_gene["mTEC","downstream"] + mTEC_gene$BarcodeCount[pos_i]
      }
    }
  }
  
  for (pos_j in 1:nrow(tissues_gene)) {
    #### keep in mind gene sense or antisense
    if (unique(strand(first_cds[i])) == '+'){
      #start and end of cds from cds objects
      start_cds <- start(first_cds[i])
      end_cds <- end(last_cds[i])
      #### counts in matrix according to position relative to start of cds
      if ((tissues_gene$start[pos_j] >= start_cds) & (tissues_gene$start[pos_j] <= end_cds) ) {
        matrix_gene["tissues","cds"] <- matrix_gene["tissues","cds"] + tissues_gene$BarcodeCount[pos_j]
      } else if(tissues_gene$start[pos_j] < start_cds) {
        matrix_gene["tissues","upstream"] <- matrix_gene["tissues","upstream"] + tissues_gene$BarcodeCount[pos_j]
      } else {
        matrix_gene["tissues","downstream"] <- matrix_gene["tissues","downstream"] + tissues_gene$BarcodeCount[pos_j]
      }

    } else {
      start_cds <- end(first_cds[i])
      end_cds <- start(last_cds[i])
      if ((tissues_gene$start[pos_j] <= start_cds) & (tissues_gene$start[pos_j] >= end_cds) )  {
        matrix_gene["tissues","cds"] <- matrix_gene["tissues","cds"] + tissues_gene$BarcodeCount[pos_j]
      } else if (tissues_gene$start[pos_j] > start_cds) {
        matrix_gene["tissues","upstream"] <- matrix_gene["tissues","upstream"] + tissues_gene$BarcodeCount[pos_j]
      } else {
        matrix_gene["tissues","downstream"] <- matrix_gene["tissues","downstream"] + tissues_gene$BarcodeCount[pos_j]
      }

    }
  }
  
  list_genes[[i]] <- matrix_gene
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

saveRDS(list_genes, file="/home/stroemic/hiwi_16/analysis/shifted_TSS/list_genes_matrices_protein_coding.rds")

