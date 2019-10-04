library(dplyr)
library(DESeq2)
library(ggplot2)

setwd("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/")
cv <- function(x) sd(x)/mean(x)

tra_entrez <- read.table(file ="/home/stroemic/hiwi_16/analysis/gene_lists/human_tra_genes_sansom.csv", sep=",", header = TRUE)
files <- list.files(pattern = "*ctss.filtered.csv$")
files <- append(files, "all_mTECs.filtered.csv")

c = 1
dats <- list()
for (i in files) {
  dat <- read.table(i, sep = ",", header=TRUE)
  if (grepl("mTEC", i)) {
    tissue <- "all_mTECs"
  } else {
    tissue <- gsub("^(.*).CNhs.*", "\\1", i )
  }
  
  print(tissue)
  #rename columns
  dat_out <- data.frame("geneID"=dat$geneID, "V2"= dat$BarcodeCount)
  #aggregate to genes 
  agg <- aggregate(dat_out$V2, by=list(geneID=dat_out$geneID), FUN=sum)
  names(agg)[names(agg) == "x"] <- tissue
  dats[[c]] <- agg
  c = c+1
}

#join all tissues together by genes
tissue_expression <- Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="geneID"), dats)
tissue_expression[is.na(tissue_expression)] <- 0
#set geneIDs as rownames and cut column away
rownames(tissue_expression) <- tissue_expression$geneID
expression_matrix <- tissue_expression[,-1]

##normalization of expression levels
sf <- estimateSizeFactorsForMatrix(expression_matrix)
expr_norm <- as.data.frame(t(t(expression_matrix)/sf))

#for_rank <- expr_norm

##calculate statistical identifiers
expr_norm$mean <- rowMeans(expr_norm)
expr_norm$median <- apply(expr_norm[,-55], 1, median)
expr_norm$sd <- apply(expr_norm[,-(55:56)], 1, sd)
expr_norm$var <- apply(expr_norm[,-(55:57)], 1, var)
expr_norm$cv <- (apply(expr_norm[,-(55:58)], 1, cv))^2

## make column with tras 
expr_norm$tra <- match(rownames(expr_norm), tra_entrez$entrezgene)
expr_norm[!is.na(expr_norm$tra),]$tra<- 'tra_in_list'
expr_norm[is.na(expr_norm$tra),]$tra <- 'not_in_list'

write.table(expr_norm,file="/home/stroemic/hiwi_16/analysis/gene_lists/tra_san_gene_expression_stats_all_tissues.csv", quote = FALSE, row.names = TRUE, col.names = TRUE, na = "", sep = ",")

pdf(file="/home/stroemic/hiwi_16/analysis/gene_lists/plots/tra_san_all_tissues_mTECs_genes_cv_over_mean.pdf")
qplot(expr_norm$mean, expr_norm$cv, color=expr_norm$tra, 
      log = "xy", size = I(0.2), xlab = "Mean", ylab = "CV²") +
      scale_color_manual(values=c("#252525", "#e31a1c")) +
      labs(color="")
dev.off()

# ####rank based
# ranks <- data.frame(t(apply(-for_rank, 1, rank, ties.method='average')))
# ranked_mTECs <- for_rank[ranks$all_mTECs == 1,]
# ranked_thymus <- for_rank[ranks$thymus == 1,]
# mTEC_specific_genes <- data.frame("entrezgene"=row.names(ranked_mTECs))
# thymus_specific_genes <- data.frame("entrezgene"=row.names(ranked_thymus))

# ##table with thymus specific genes 
# write.table(thymus_specific_genes,file ="/home/stroemic/hiwi_16/analysis/gene_lists/thymus_specific_genes_by_rank.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
# ##table with mTEC specific genes
# write.table(mTEC_specific_genes,file ="/home/stroemic/hiwi_16/analysis/gene_lists/mTEC_specific_genes_by_rank.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")


# ####glmnet
# for_glmnet$is_mTEC <- "blub"
# for_glmnet[(for_glmnet$all_mTECs == 0),]$is_mTEC<- 'no'
# for_glmnet[!(for_glmnet$all_mTECs == 0),]$is_mTEC<- 'yes'
# 
# test <- glmnet(as.matrix(for_glmnet[,-(54:55)]), for_glmnet$is_mTEC, family = "binomial")
# cv <- cv.glmnet(as.matrix(for_glmnet[,-(54:55)]), for_glmnet$is_mTEC, family = "binomial")
# 
# 
# ####DESeq2
# tissue_expression[,-1]<- sapply(round(tissue_expression[,-1], 0), as.integer)
# counts <- as.matrix(tissue_expression[,-1])
# columns <- data.frame(row.names = colnames(tissue_expression[,-1]), "tissue" = colnames(tissue_expression[,-1]), "dataset"=c("Fantom5","mTEC"))
# columns[!(columns$tissue == "all_mTECs"),]$dataset <- "Fantom5"
# columns <- subset(columns, select=-c(tissue))
# dds <- DESeqDataSetFromMatrix(countData = counts, colData = columns, design = ~ dataset)
# start_time <- Sys.time()
# dds <- DESeq(dds, parallel = TRUE)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# res <- results(dds, contrast = c("dataset","Fantom5","mTEC"))
# filtered_res <- res[res$log2FoldChange < (-3.00) & res$pvalue < 0.05,]
# 

