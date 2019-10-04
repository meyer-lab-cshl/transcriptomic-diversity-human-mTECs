library(dplyr)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)
library(reshape2)
library(TTR)

#### functions ####
preprocess_df <- function(x) {
  x <- x[x$geneID %in% names(tx_genes),]
  x <- droplevels(x)
  x <- x[with(x, order(chrom, geneID, start)),]
  filt <- x %>% group_by(geneID) %>% mutate(the_rank  = rank(-BarcodeCount, ties.method = "first")) %>% filter(the_rank == 1)
  
  x$geneID <- as.character(x$geneID)
  
  x_plus <- x[x$strand == "+",]
  x_minus <- x[x$strand == "-",]
  
  #x_plus$rel_pos <- with(x_plus, (x_plus$start - filt$start[match(x_plus$geneID, filt$geneID)]))
  #x_minus$rel_pos <- with(x_minus, ( filt$start[match(x_minus$geneID, filt$geneID)] - x_minus$start))
  
  x_plus$rel_pos <- with(x_plus, (x_plus$start - start(tx_genes[x_plus$geneID])))
  x_minus$rel_pos <- with(x_minus, ( end(tx_genes[x_minus$geneID]) - x_minus$start))
  
  x <- rbind(x_plus, x_minus)
  x$rel_pos <- as.character(x$rel_pos)
  return(x)
}

plot_range <- function(m_df) {
  plot <- ggplot(m_df, aes(x=Var2, y=value, group =Var1, colour=Var1 )) + 
                  geom_point(size=0.5) +
                  theme_bw() +
                  scale_x_continuous(limits = c(-40, 40)) +
                  labs(title="Distribution of tags over genomic range for different gene classes", x ="bp around annotated TSS", y ="", colour="")
  return(plot)
}

plot_range_line <- function(m_df) {
  plot <- ggplot(m_df, aes(x=Var2, y=value, group =Var1, colour=Var1 )) + 
                  geom_point(size=0.5) +
                  theme_bw() +
                  geom_line(aes(y=roll)) +
                  scale_x_continuous(limits = c(-40, 40)) +
                  #geom_smooth(method="loess", se=F, span=0.05) +
                  labs(title="Distribution of tags over genomic range for different gene classes", x ="bp around annotated TSS", y ="", colour="")
  return(plot)
}

plot_range_line_only <- function(m_df) {
  plot <- ggplot(m_df, aes(x=Var2, y=value, group =Var1, colour=Var1 )) + 
                  #geom_point(size=0.5) +
                  theme_bw() +
                  geom_line(aes(y=roll)) +
                  scale_x_continuous(limits = c(-40, 40)) +
                  #geom_smooth(method="loess", se=F, span=0.05) +
                  labs(title="Distribution of tags over genomic range for different gene classes", x ="bp around annotated TSS", y ="", colour="")
  return(plot)
}

#### read in of data frames ####
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
mTEC_pos_fil <- read.table("/home/stroemic/hiwi_16/data/raw_positions/all_mTECs.positions.csv", header=TRUE, sep=",")
##tissue (wo thymus) all positions
tissues_pos_fil <- read.table("/home/stroemic/hiwi_16/data/raw_positions/all.tissues.wo.thymus.positions.csv", header = TRUE, sep = ",")


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


#### Analysis ####
genes_entrez_mTECs <- list(aire_dep_tras$entrezgene, fezf2_dep_tras$entrezgene, other_tras$entrezgene, hk_genes$entrezgene ,mTEC_genes$entrezgene, thymus_genes$entrezgene)
genes_entrez_tissues <- list(hk_genes$entrezgene, aire_dep_tras$entrezgene, fezf2_dep_tras$entrezgene, other_tras$entrezgene)

gene_classes_mTECs <- c("aire_dep_tras_mTECs", "fezf2_dep_tras_mTECs", "other_tras_mTECs", "hk_mTECs", "mTEC_spec_mTECs", "thymus_spec_mTECs")
gene_classes_tissues <- c("hk_tissues", "aire_dep_tras_tissues", "fezf2_dep_tras_tissues", "other_tras_tissues")

mTEC_pos_fil <- preprocess_df(mTEC_pos_fil)
tissues_pos_fil <- preprocess_df(tissues_pos_fil)

list_uncounted_pos <- list()
list_matrices <- list()
list_counts_classes <- list()

for (i in c(250, 1000, 20000) ) {
  name <- as.character(i)
  if (i == 250) {
    thres=-250
    counts_mTECs <- matrix(0, nrow = length(unique(mTEC_pos_fil$geneID)), ncol = (i+251), 
                           dimnames = list(unique(mTEC_pos_fil$geneID), as.character(c(-250:i))))
    
    counts_tissues <- matrix(0, nrow = length(unique(tissues_pos_fil$geneID)), ncol = (i+251),
                             dimnames = list(unique(tissues_pos_fil$geneID), as.character(c(-250:i))))
  } else {
    thres= -1000
    counts_mTECs <- matrix(0, nrow = length(unique(mTEC_pos_fil$geneID)), ncol = (i+1001), 
                           dimnames = list(unique(mTEC_pos_fil$geneID), as.character(c(-1000:i))))
    
    counts_tissues <- matrix(0, nrow = length(unique(tissues_pos_fil$geneID)), ncol = (i+1001),
                             dimnames = list(unique(tissues_pos_fil$geneID), as.character(c(-1000:i))))
  }

  uncounted_pos_m = 0
  start.time <- Sys.time()
  print("calculating mTEC matrix")
  for (x in 1:nrow(mTEC_pos_fil)) {
    rel_position <- mTEC_pos_fil[x,]$rel_pos
    gene <- mTEC_pos_fil[x,]$geneID
    if (as.integer(rel_position) > i | as.integer(rel_position) < thres) {
      uncounted_pos_m = uncounted_pos_m + 1
    } else  {
      counts_mTECs[gene, rel_position] <- (counts_mTECs[gene, rel_position] + mTEC_pos_fil[x,]$BarcodeCount)
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  
  uncounted_pos_t = 0
  start.time <- Sys.time()
  for (y in 1:nrow(tissues_pos_fil)) {
    rel_position <- tissues_pos_fil[y,]$rel_pos
    gene <- tissues_pos_fil[y,]$geneID
    if (as.integer(rel_position) > i | as.integer(rel_position) < thres) {
      uncounted_pos_t = uncounted_pos_t + 1
    } else  {
      counts_tissues[gene, rel_position] <- (counts_tissues[gene, rel_position] + tissues_pos_fil[y,]$BarcodeCount)
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  list_uncounted_pos[[name]] <-list("uncounted_pos_m"=uncounted_pos_m, "uncounted_pos_t"=uncounted_pos_t) 
  list_matrices[[name]] <-list("counts_mTECs"=counts_mTECs, "counts_tissues"=counts_tissues)
  
  list_mTECs <- list()
  
  for (j in 1:6) {
    print(gene_classes_mTECs[j])
    counts_j <- subset(list_matrices[[name]][[1]], rownames(list_matrices[[name]][[1]]) %in% genes_entrez_mTECs[[j]])
    print(dim(counts_j))
    sums <- colSums(counts_j)
    sums_norm <- (sums-min(sums))/(max(sums)-min(sums))
    list_mTECs[[gene_classes_mTECs[j]]] <- sums_norm
  }
  
  list_tissues <- list()
  
  for (k in 1:4) {
    print(gene_classes_tissues[k])
    counts_k <- subset(list_matrices[[name]][[2]], rownames(list_matrices[[name]][[2]]) %in% genes_entrez_tissues[[k]])
    print(dim(counts_k))
    sums <- colSums(counts_k)
    sums_norm <- (sums-min(sums))/(max(sums)-min(sums))
    list_tissues[[gene_classes_tissues[k]]] <- sums_norm
  }
  
  list_both <- c(list_mTECs, list_tissues)
  counts_classes <- do.call(rbind, list_both)
  list_counts_classes[[name]] <- counts_classes
}

saveRDS(list_uncounted_pos, file ="/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_rel_pos_to_annotation_list_uncounted_positions_human_filt.rds")
saveRDS(list_matrices, file="/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_rel_pos_to_annotation_list_counts_matrices_human_filt.rds")
saveRDS(list_counts_classes, file="/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_rel_pos_to_annotation_list_counts_classes_human_filt.rds")

# saveRDS(list_uncounted_pos, file =paste("/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_rel_pos_to_annotation_list_uncounted_positions_human_filt_",name,".rds",sep=""))
# saveRDS(list_matrices, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_rel_pos_to_annotation_list_counts_matrices_human_filt_",name,".rds",sep=""))
# saveRDS(list_counts_classes, file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_rel_pos_to_annotation_list_counts_classes_human_filt_",name,".rds",sep=""))

# ###when loaded
#  counts_250 <- readRDS(file="/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_list_counts_classes_human_filt_250.rds")
#  counts_1000 <- readRDS(file="/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_list_counts_classes_human_filt_1000.rds")
#  counts_20000 <- readRDS(file="/home/stroemic/hiwi_16/analysis/gene_lists/raw_data_list_counts_classes_human_filt_20000.rds")
#  list_counts_classes <- c(counts_250, counts_1000, counts_20000)
# 
# 
# #### plotting ####
# for (i in c(250)) {
#   name <- as.character(i)
#   m_counts_classes_tra_m <- melt(list_counts_classes[[name]][c(1,2,3),])
#   m_counts_classes_tra_m <- m_counts_classes_tra_m %>% group_by(Var1) %>% mutate("roll" =SMA(value, n=5))
#   m_counts_classes_aire_tra <- melt(list_counts_classes[[name]][c(1,8),])
#   m_counts_classes_aire_tra <- m_counts_classes_aire_tra %>% group_by(Var1) %>% mutate("roll" =SMA(value, n=5))
#   m_counts_classes_fezf2_tra <- melt(list_counts_classes[[name]][c(2,9),])
#   m_counts_classes_fezf2_tra <- m_counts_classes_fezf2_tra %>% group_by(Var1) %>% mutate("roll" =SMA(value, n=5))
#   m_counts_classes_other_tra <- melt(list_counts_classes[[name]][c(3,10),])
#   m_counts_classes_other_tra <- m_counts_classes_other_tra%>% group_by(Var1) %>% mutate("roll" =SMA(value, n=5))
#   m_counts_classes_hk <- melt(list_counts_classes[[name]][c(4,7),])
#   m_counts_classes_hk <- m_counts_classes_hk %>% group_by(Var1) %>% mutate("roll" =SMA(value, n=5))
#   m_counts_classes_tra_hk_m <- melt(list_counts_classes[[name]][c(1:4),])
#   m_counts_classes_tra_hk_m <- m_counts_classes_tra_hk_m%>% group_by(Var1) %>% mutate("roll" =SMA(value, n=5))
#   m_counts_classes_m_t_m <- melt(list_counts_classes[[name]][c(5,6),])
#   m_counts_classes_m_t_m <- m_counts_classes_m_t_m%>% group_by(Var1) %>% mutate("roll" =SMA(value, n=5))
#   
# 
#   pdf(file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/plots/raw_data_gene_range_analysis_",name,"_40.pdf",sep=""))
#   print(plot <- plot_range(m_counts_classes_tra_m))
#   print(plot <- plot_range(m_counts_classes_aire_tra))
#   print(plot <- plot_range(m_counts_classes_fezf2_tra))
#   print(plot <- plot_range(m_counts_classes_other_tra))
#   print(plot <- plot_range(m_counts_classes_hk))
#   print(plot <- plot_range(m_counts_classes_tra_hk_m))
#   print(plot <- plot_range(m_counts_classes_m_t_m))
#   dev.off()
#   pdf(file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/plots/raw_data_gene_range_analysis_line_",name,"_40.pdf",sep=""))
#   print(plot <- plot_range_line(m_counts_classes_tra_m))
#   print(plot <- plot_range_line(m_counts_classes_aire_tra))
#   print(plot <- plot_range_line(m_counts_classes_fezf2_tra))
#   print(plot <- plot_range_line(m_counts_classes_other_tra))
#   print(plot <- plot_range_line(m_counts_classes_hk))
#   print(plot <- plot_range_line(m_counts_classes_tra_hk_m))
#   print(plot <- plot_range_line(m_counts_classes_m_t_m))
#   dev.off()
#   pdf(file=paste("/home/stroemic/hiwi_16/analysis/gene_lists/plots/raw_data_gene_range_analysis_line_only_",name,"_40.pdf",sep=""))
#   print(plot <- plot_range_line_only(m_counts_classes_tra_m))
#   print(plot <- plot_range_line_only(m_counts_classes_aire_tra))
#   print(plot <- plot_range_line_only(m_counts_classes_fezf2_tra))
#   print(plot <- plot_range_line_only(m_counts_classes_other_tra))
#   print(plot <- plot_range_line_only(m_counts_classes_hk))
#   print(plot <- plot_range_line_only(m_counts_classes_tra_hk_m))
#   print(plot <- plot_range_line_only(m_counts_classes_m_t_m))
#   dev.off()
# }
#  
