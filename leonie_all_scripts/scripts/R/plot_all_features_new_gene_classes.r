library(dplyr)
library(raster)
library(ggplot2)
library(reshape2)
library(biomaRt)
library(scales)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

#### read in of data frames ####
##table with tras
tra_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_tra_genes_sansom.csv", header=TRUE, sep=",")
##table with aire dep genes 
aire_dep_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_aire_dep_san_genes.csv", header=TRUE, sep=",")
##table with fezf2 dep genes 
fezf2_dep_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_fezf2_dep_genes.csv", header=TRUE, sep=",")
##table with housekeeping genes defined on all tissues
hk_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/housekeeping_genes.csv", header=TRUE, sep=",")

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

# ##Filter input data to contain only protein coding genes
# mTEC_pos_fil <- mTEC_pos_fil[mTEC_pos_fil$geneID %in% names(tx_genes),]
# mTEC_pos_fil <- droplevels(mTEC_pos_fil)
# tissues_pos_fil <- tissues_pos_fil[tissues_pos_fil$geneID %in% names(tx_genes),]
# tissues_pos_fil <- droplevels(tissues_pos_fil)
# mTEC_pos_unfil <- mTEC_pos_unfil[mTEC_pos_unfil$geneID %in% names(tx_genes),]
# mTEC_pos_unfil <- droplevels(mTEC_pos_unfil)
# tissues_pos_unfil <- tissues_pos_unfil[tissues_pos_unfil$geneID %in% names(tx_genes),]
# tissues_pos_unfil <- droplevels(tissues_pos_unfil)

features_dfs <- list(human_filt=list(mTEC_pos_fil, tissues_pos_fil), 
                     human_unfilt=list(mTEC_pos_unfil, tissues_pos_unfil))

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

#### Define GRanges for Overlaps ####
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
# saveRDS(introns, file="/home/stroemic/hiwi_16/analysis/gene_lists/introns.rds")
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

for (j in 1:2) {
  class <- names(features_dfs)[j]
  print(class)
  #### Analysis for mTEC data ####
  ## make range objects for tags
  mTEC_ranges <- makeGRangesFromDataFrame(features_dfs[[j]][[1]], keep.extra.columns = TRUE)
  
  gene_classes_mTECs <- c("aire_dep_tras_mTECs", "fezf2_dep_tras_mTECs", "other_tras_mTECs", "hk_mTECs", "other_genes_mTECs")
  genes_entrez_mTECs <- list(aire_dep_tras$entrezgene, fezf2_dep_tras$entrezgene, other_tras$entrezgene, hk_genes$entrezgene, others )
  
  
  ## findOverlaps and store in vectors
  TSS10_overlaps <- findOverlaps(mTEC_ranges, TSS10_genes, select = "first")
  TSS100_overlaps <- findOverlaps(mTEC_ranges, TSS100_genes, select = "first")
  first_exon_overlaps <- findOverlaps(mTEC_ranges, first_exon, select = "first")
  other_exon_overlaps <- findOverlaps(mTEC_ranges, other_exon, select = "first")
  TSS10_oe_overlaps <- findOverlaps(mTEC_ranges, TSS10_other_exon, select = "first")
  TSS100_oe_overlaps <- findOverlaps(mTEC_ranges, TSS100_other_exon, select = "first")
  intron_overlaps <- findOverlaps(mTEC_ranges, introns, select = "first")
  TTS_100_overlaps <- findOverlaps(mTEC_ranges, TTS100_genes, select = "first")
  upstream_overlaps <- findOverlaps(mTEC_ranges, upstream_genes, select = "first")
  downstream_overlaps <- findOverlaps(mTEC_ranges, downstream_genes, select = "first")
  antisense_overlaps <- findOverlaps(mTEC_ranges, antisense, select = "first")
  
  ## define matrix to store counts calculated in loop
  counts_mTECs <- matrix(0, nrow=12, ncol=6) 
  rownames(counts_mTECs) <-c("TSS_10","TSS_100","first_exon","TSS_10_other_exons","TSS_100_other_exons","TTS_100","upstream","other_exon","intron","antisense","downstream","intergenic") 
  colnames(counts_mTECs) <-c("aire_dep_tras_mTECs", "fezf2_dep_tras_mTECs", "other_tras_mTECs", "hk_mTECs", "other_genes_mTECs","rest_mTECs")                
  
  start.time <- Sys.time()
  print("mTEC matrix")
  for (i in 1:length(mTEC_ranges) ) { 
    if (!is.na(TSS10_overlaps[i])){
      feature_class <- "TSS_10"
      gene <- names(TSS10_genes[TSS10_overlaps[i]])
      
    } else if(!is.na(TSS100_overlaps[i])) {
      feature_class <- "TSS_100"
      gene <- names(TSS100_genes[TSS100_overlaps[i]])
      
    } else if(!is.na(first_exon_overlaps[i])) {
      feature_class <- "first_exon"
      gene <- names(first_exon[first_exon_overlaps[i]])
      
    } else if(!is.na(TSS10_oe_overlaps[i])) {
      feature_class <- "TSS_10_other_exons"
      gene <- names(TSS10_other_exon[TSS10_oe_overlaps[i]])
      
    } else if(!is.na(TSS100_oe_overlaps[i])) {
      feature_class <- "TSS_100_other_exons"
      gene <- names(TSS100_other_exon[TSS100_oe_overlaps[i]])
      
    } else if(!is.na(TTS_100_overlaps[i])) {
      feature_class <- "TTS_100"
      gene <- names(TTS100_genes[TTS_100_overlaps[i]])
      
    } else if(!is.na(upstream_overlaps[i])) {
      feature_class <- "upstream"
      gene <- names(upstream_genes[upstream_overlaps[i]])
      
    } else if(!is.na(other_exon_overlaps[i])) {
      feature_class <- "other_exon"
      gene <- names(other_exon[other_exon_overlaps[i]])
      
    } else if(!is.na(intron_overlaps[i])) {
      feature_class <- "intron"
      gene <- names(introns[intron_overlaps[i]])
      
    } else if(!is.na(antisense_overlaps[i])) {
      feature_class <- "antisense"
      gene <- names(antisense[antisense_overlaps[i]])
      
    } else if(!is.na(downstream_overlaps[i])) {
      feature_class <- "downstream"
      gene <- names(downstream_genes[downstream_overlaps[i]])
      
    } else { 
      feature_class <- "intergenic" 
      gene <- NA
    }
    
    
    
    l <- sapply(genes_entrez_mTECs, function(x) is.element(gene,x) )
    in_gene_classes <- gene_classes_mTECs[l] 
    if (identical(in_gene_classes, character(0))) in_gene_classes <- "rest_mTECs"
    
    counts_mTECs[feature_class, in_gene_classes] <- (counts_mTECs[feature_class, in_gene_classes] + features_dfs[[j]][[1]]$BarcodeCount[i])
    
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  saveRDS(counts_mTECs, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/new_classes_tss_oe_tra_san_aire_san_counts_mTECs_matrix_",class,".rds", sep=""))
  print("saved")

  #### Analysis for tissue data ####
  ## make range objects for tags
  tissue_ranges <- makeGRangesFromDataFrame(features_dfs[[j]][[2]], keep.extra.columns = TRUE)
  
  gene_classes_tissues <- c("aire_dep_tras_tissues", "fezf2_dep_tras_tissues", "other_tras_tissues", "hk_tissues","other_genes_tissues")    ####change here
  genes_entrez_tissues <- list(aire_dep_tras$entrezgene, fezf2_dep_tras$entrezgene, other_tras$entrezgene, hk_genes$entrezgene, others)    ####change here
  
  
  
  ## findOverlaps and store in vectors
  TSS10_overlaps <- findOverlaps(tissue_ranges, TSS10_genes, select = "first")
  TSS100_overlaps <- findOverlaps(tissue_ranges, TSS100_genes, select = "first")
  first_exon_overlaps <- findOverlaps(tissue_ranges, first_exon, select = "first")
  other_exon_overlaps <- findOverlaps(tissue_ranges, other_exon, select = "first")
  TSS10_oe_overlaps <- findOverlaps(tissue_ranges, TSS10_other_exon, select = "first")
  TSS100_oe_overlaps <- findOverlaps(tissue_ranges, TSS100_other_exon, select = "first")
  intron_overlaps <- findOverlaps(tissue_ranges, introns, select = "first")
  TTS_100_overlaps <- findOverlaps(tissue_ranges, TTS100_genes, select = "first")
  upstream_overlaps <- findOverlaps(tissue_ranges, upstream_genes, select = "first")
  downstream_overlaps <- findOverlaps(tissue_ranges, downstream_genes, select = "first")
  antisense_overlaps <- findOverlaps(tissue_ranges, antisense, select = "first")
  
  ## define matrix to store counts calculated in loop
  counts_tissues <- matrix(0, nrow=12, ncol=6) ####change here
  rownames(counts_tissues) <-c("TSS_10","TSS_100","first_exon","TSS_10_other_exons","TSS_100_other_exons","TTS_100","upstream","other_exon","intron","antisense","downstream","intergenic") ##intron
  colnames(counts_tissues) <-c( "aire_dep_tras_tissues", "fezf2_dep_tras_tissues", "other_tras_tissues", "hk_tissues", "other_genes_tissues","rest_tissues")   ####change here
  
  
  start.time <- Sys.time()
  print("tissues matrix")
  for (i in 1:length(tissue_ranges)) { 
    if (!is.na(TSS10_overlaps[i])){
      feature_class <- "TSS_10"
      gene <- names(TSS10_genes[TSS10_overlaps[i]])
      
    } else if(!is.na(TSS100_overlaps[i])) {
      feature_class <- "TSS_100"
      gene <- names(TSS100_genes[TSS100_overlaps[i]])
      
    } else if(!is.na(first_exon_overlaps[i])) {
      feature_class <- "first_exon"
      gene <- names(first_exon[first_exon_overlaps[i]])
      
    } else if(!is.na(TSS10_oe_overlaps[i])) {
      feature_class <- "TSS_10_other_exons"
      gene <- names(TSS10_other_exon[TSS10_oe_overlaps[i]])
      
    } else if(!is.na(TSS100_oe_overlaps[i])) {
      feature_class <- "TSS_100_other_exons"
      gene <- names(TSS100_other_exon[TSS100_oe_overlaps[i]])
      
    } else if(!is.na(TTS_100_overlaps[i])) {
      feature_class <- "TTS_100"
      gene <- names(TTS100_genes[TTS_100_overlaps[i]])
      
    } else if(!is.na(upstream_overlaps[i])) {
      feature_class <- "upstream"
      gene <- names(upstream_genes[upstream_overlaps[i]])
      
    } else if(!is.na(other_exon_overlaps[i])) {
      feature_class <- "other_exon"
      gene <- names(other_exon[other_exon_overlaps[i]])
      
    } else if(!is.na(intron_overlaps[i])) {
      feature_class <- "intron"
      gene <- names(introns[intron_overlaps[i]])
      
    } else if(!is.na(antisense_overlaps[i])) {
      feature_class <- "antisense"
      gene <- names(antisense[antisense_overlaps[i]])
      
    } else if(!is.na(downstream_overlaps[i])) {
      feature_class <- "downstream"
      gene <- names(downstream_genes[downstream_overlaps[i]])
      
    } else { 
      feature_class <- "intergenic" 
      gene <- NA
    }
    
    
    
    l <- sapply(genes_entrez_tissues, function(x) is.element(gene,x) )
    in_gene_classes <- gene_classes_tissues[l] 
    if (identical(in_gene_classes, character(0))) in_gene_classes <- "rest_tissues"
    
    counts_tissues[feature_class, in_gene_classes] <- (counts_tissues[feature_class, in_gene_classes] + features_dfs[[j]][[2]]$BarcodeCount[i])
    
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  saveRDS(counts_tissues, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/new_classes_tss_oe_tra_san_aire_san_counts_tissues_matrix_",class,".rds", sep=""))
  print("saved")

  #### further calculations and plotting ####
  counts <- cbind(counts_mTECs, counts_tissues)
  counts_trunc <- counts[-12,] ####change here
  m_counts <- melt(counts_trunc)
  m_counts <- transform(m_counts, feature_class = factor(Var1, rownames(counts_trunc)),
                        gene_class = factor(Var2, colnames(counts_trunc)))
  
  pdf(file =paste("/home/stroemic/hiwi_16/analysis/gene_lists/plots/newc_classes_tss_oe_tra_san_aire_san_absolut_positions_proportions_by_matrix_",class,".pdf", sep=""))
  print(ggplot(m_counts, aes(x=gene_class,y=value,fill=feature_class)) + 
          geom_bar(position = "fill", stat = "identity") +
          labs(title="Porportion of tags belonging to each feature in different gene sets", x ="", y ="Proportion") +
          theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
          scale_fill_manual(values = c('#a6cee3','#1f78b4','#fb9a99','#b2df8a','#33a02c','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99')) +
          guides(fill=guide_legend(title=NULL))
        
  )
  dev.off()
  
  
  pdf(file =paste("/home/stroemic/hiwi_16/analysis/gene_lists/plots/new_classes_tss_oe_tra_san_aire_san_absolut_positions_percentage_by_matrix_",class,".pdf",sep=""))
  counts_trunc <- t(t(counts_trunc)/rowSums(t(counts_trunc)))
  for (i in 1:11) {
    name <- rownames(counts_trunc)[i]
    print(name)
    ##define df for plotting using percentage function
    plot_df <- data.frame("feature"=counts_trunc[i,], 
                          "source"=factor(names(counts_trunc[i,]), levels = c("aire_dep_tras_mTECs", "fezf2_dep_tras_mTECs", "other_tras_mTECs", "hk_mTECs","other_genes_mTECs","rest_mTECs", "aire_dep_tras_tissues", "fezf2_dep_tras_tissues", "other_tras_tissues", "hk_tissues", "other_genes_tissues","rest_tissues")))
    print(ggplot(plot_df, aes(x= source, y= feature)) + 
            geom_bar(stat = "identity", fill="#1f78b4") +
            labs(title=paste("Percentage of tags belonging to feature \"",name,"\" in different gene sets", sep = ""), x ="", y =name ) +
            theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
            scale_y_continuous(labels = percent)
    )
  }
  dev.off()
}
