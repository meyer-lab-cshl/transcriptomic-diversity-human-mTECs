library(dplyr)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(biomaRt)
library(reshape2)

#### functions ####
preprocess_df <- function(x, keep) {
  x <- x[x$geneID %in% names(first_cds),]
  x <- x[x$geneID %in% keep,]
  x <- droplevels(x)
  x$geneID <- as.character(x$geneID)
  return(x)
}

#### read in ####
##internal all positions
internal_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_internal_mm9_ESC.filtered.csv", header=TRUE, sep=",")
##tissue (wo thymus) all positions
fantom5_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.fantom5.mm9.ESC.filtered.csv", header = TRUE, sep = ",")
##get common genes for calculation
keep <- intersect(internal_pos_fil$geneID, fantom5_pos_fil$geneID)

##### BioMart object for mouse version NCBIM37 = mm9
ensembl_mm9=useMart(host='may2012.archive.ensembl.org', 
                  biomart='ENSEMBL_MART_ENSEMBL', 
                  dataset="mmusculus_gene_ensembl")

##### txdb gene object
my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19",  "chrM", "chrX", "chrY")
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)
### keep only protein coding
mart_gene_type <- getBM(attributes = c('entrezgene', 'gene_biotype'), 
                        filters = 'entrezgene',
                        values = names(tx_genes), 
                        mart=ensembl_mm9)
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

internal_pos_fil <- preprocess_df(internal_pos_fil, keep)
fantom5_pos_fil <- preprocess_df(fantom5_pos_fil, keep)

list_genes_names <- names(first_cds)
list_genes <- vector("list", length(list_genes_names))
names(list_genes) <- list_genes_names

start.time <- Sys.time()
c=1
for (i in names(list_genes)) {
  print(c)
  c=c+1
  #build matrix to write counts into 
  matrix_gene <- matrix(0, nrow=2, ncol=3)
  rownames(matrix_gene) <-c("fantom5","internal") 
  colnames(matrix_gene) <-c("upstream","cds","downstream")
  #subdata frames for the gene
  internal_gene <- internal_pos_fil[internal_pos_fil$geneID == i,]
  fantom5_gene <- fantom5_pos_fil[fantom5_pos_fil$geneID == i,]
  #iterate over all positions in internal data for the gene 
  for (pos_i in 1:nrow(internal_gene)) {
    #### keep in mind gene sense or antisense
    if (unique(strand(first_cds[i])) == '+'){
      #start and end of cds from cds objects
      start_cds <- start(first_cds[i])
      end_cds <- end(last_cds[i])
      #### counts in matrix according to position either in cds, upstream or downstream
      if ((internal_gene$start[pos_i] >= start_cds) & (internal_gene$start[pos_i] <= end_cds)) {
        matrix_gene["internal","cds"] <- matrix_gene["internal","cds"] + internal_gene$BarcodeCount[pos_i]
      } else if (internal_gene$start[pos_i] < start_cds) {
        matrix_gene["internal","upstream"] <- matrix_gene["internal","upstream"] + internal_gene$BarcodeCount[pos_i]
      } else {
        matrix_gene["internal","downstream"] <- matrix_gene["internal","downstream"] + internal_gene$BarcodeCount[pos_i]
      }
      
    } else {
      ## negative strand with different start and end 
      start_cds <- end(first_cds[i])
      end_cds <- start(last_cds[i])
      if ((internal_gene$start[pos_i] <= start_cds) & (internal_gene$start[pos_i] >= end_cds)) {
        matrix_gene["internal","cds"] <- matrix_gene["internal","cds"] + internal_gene$BarcodeCount[pos_i]
      } else if (internal_gene$start[pos_i] > start_cds) {
        matrix_gene["internal","upstream"] <- matrix_gene["internal","upstream"] + internal_gene$BarcodeCount[pos_i]
      } else { 
        matrix_gene["internal","downstream"] <- matrix_gene["internal","downstream"] + internal_gene$BarcodeCount[pos_i]
      }
    }
  }
  
  for (pos_j in 1:nrow(fantom5_gene)) {
    #### keep in mind gene sense or antisense
    if (unique(strand(first_cds[i])) == '+'){
      #start and end of cds from cds objects
      start_cds <- start(first_cds[i])
      end_cds <- end(last_cds[i])
      #### counts in matrix according to position relative to start of cds
      if ((fantom5_gene$start[pos_j] >= start_cds) & (fantom5_gene$start[pos_j] <= end_cds) ) {
        matrix_gene["fantom5","cds"] <- matrix_gene["fantom5","cds"] + fantom5_gene$BarcodeCount[pos_j]
      } else if(fantom5_gene$start[pos_j] < start_cds) {
        matrix_gene["fantom5","upstream"] <- matrix_gene["fantom5","upstream"] + fantom5_gene$BarcodeCount[pos_j]
      } else {
        matrix_gene["fantom5","downstream"] <- matrix_gene["fantom5","downstream"] + fantom5_gene$BarcodeCount[pos_j]
      }
      
    } else {
      start_cds <- end(first_cds[i])
      end_cds <- start(last_cds[i])
      if ((fantom5_gene$start[pos_j] <= start_cds) & (fantom5_gene$start[pos_j] >= end_cds) )  {
        matrix_gene["fantom5","cds"] <- matrix_gene["fantom5","cds"] + fantom5_gene$BarcodeCount[pos_j]
      } else if (fantom5_gene$start[pos_j] > start_cds) {
        matrix_gene["fantom5","upstream"] <- matrix_gene["fantom5","upstream"] + fantom5_gene$BarcodeCount[pos_j]
      } else {
        matrix_gene["fantom5","downstream"] <- matrix_gene["fantom5","downstream"] + fantom5_gene$BarcodeCount[pos_j]
      }
      
    }
  }
  
  list_genes[[i]] <- matrix_gene
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

saveRDS(list_genes, file="/home/stroemic/hiwi_16/analysis/shifted_TSS/list_genes_matrices_mm9_protein_coding.rds")
# 
# list_genes_mm9 <- readRDS(file="/home/stroemic/hiwi_16/analysis/shifted_TSS/list_genes_matrices_mm9.rds")
# list_fisher <- lapply(list_genes_mm9, function(x) fisher.test(round(x[,1:2])))
# ratio_fantom5 <- sapply(list_genes_mm9, function(x) x[1,2] / (x[1,1] + x[1,2]))
# ratio_internal  <- sapply(list_genes_mm9, function(x) x[2,2] / (x[2,1] + x[2,2]))
# 
# df_ratios <- data.frame("fantom5"=ratio_fantom5, "internal"=ratio_internal)
# m_df_ratios <- melt(df_ratios)
# 
# pdf(file="/home/stroemic/hiwi_16/analysis/shifted_TSS/plots/boxplot_ratios_mm9_fantom5_internal.pdf")
# p1 <- ggplot(m_df_ratios, aes(x=variable, y=value, fill=variable))
# ylim1 = boxplot.stats(m_df_ratios$value)$stats[c(1, 5)]
# p1 + geom_boxplot(outlier.size=NA) +  coord_cartesian(ylim = ylim1*1.8)
# dev.off()
