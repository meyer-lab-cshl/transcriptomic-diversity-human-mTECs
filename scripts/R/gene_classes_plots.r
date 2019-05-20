library(dplyr)
library(raster)
library(ggplot2)
library(reshape2)
library(biomaRt)
library(scales)

percentage <- function(x,i) sum(x[x[[i]]==TRUE,]$BarcodeCount)/sum(x$BarcodeCount)
prepare_df <- function(x,i,list) {
  x[list][x[i] == TRUE,] <- FALSE
  return(x)
}

####
##table with tras
tra_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_tra_genes.csv", header=TRUE, sep=",")
##table with aire dep genes 
aire_dep_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_aire_dep_genes.csv", header=TRUE, sep=",")
##table with fezf2 dep genes 
fezf2_dep_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_fezf2_dep_genes.csv", header=TRUE, sep=",")
##table with thymus specific genes 
thymus_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/thymus_specific_genes.csv", header=TRUE, sep=",")
##table with mTEC specific genes
mTEC_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/mTEC_specific_genes.csv", header=TRUE, sep=",")
##table with housekeeping genes defined on all tissues
hk_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/housekeeping_genes.csv", header=TRUE, sep=",")

##### human 
###unfiltered
##mtec features 
mTEC_features <- read.table("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/features/all_mTECs.features.csv", header=TRUE, sep=",")
##tissue (wo thymus) features
tissues_features <- read.table("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/features/all.tissues.wo.thymus.features.csv", header = TRUE, sep = ",")

###filtered
##mtec features 
mTEC_features_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.features-filtered.csv", header=TRUE, sep=",")
mTEC_features_fil <- mTEC_features_fil[,-33]
##tissue (wo thymus) features
tissues_features_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.wo.thymus.features-filtered.csv", header = TRUE, sep = ",")
tissues_features_fil <- tissues_features_fil[,-33]


#### prepare feature tables to get rid of overlaps and counts add up to 1
dfs <- list(mTEC_features, tissues_features, mTEC_features_fil, tissues_features_fil)
dfs_prepared <- list()
c <- 1
for (df in dfs) {
  column_ids <- c(10,11,14,17,15)
  ids <- c(10,11,14,17,15,16)
  for (i in column_ids){
    ids <- ids[!ids %in% i]
    df <- prepare_df(df, i, ids)
  }
  dfs_prepared[[c]] <- df
  c = c+1
}


features_dfs <- list(human_unfilt=list(dfs_prepared[[1]], dfs_prepared[[2]]), 
                     human_filt=list(dfs_prepared[[3]], dfs_prepared[[4]]))

for (i in 1:2) {
  class <- names(features_dfs)[i]
  print(class)
  ###define datatables
  # aire dep tras in mTECs
  aire_dep_tras <- merge(tra_genes, aire_dep_genes, by="entrezgene")
  aire_dep_tras_mTECs <- merge(aire_dep_tras, features_dfs[[i]][[1]], by.x="entrezgene", by.y="geneID")
  aire_dep_tras_mTECs <- aire_dep_tras_mTECs[,-(2:3)]
  aire_dep_tras_mTECs$source <- "aire_dep_tras_mTECs"
  # fezf2 dep tras in mTECs
  fezf2_dep_tras <- merge(tra_genes, fezf2_dep_genes, by="entrezgene")
  fezf2_dep_tras_mTECs <- merge(fezf2_dep_tras, features_dfs[[i]][[1]], by.x="entrezgene", by.y="geneID")
  fezf2_dep_tras_mTECs <- fezf2_dep_tras_mTECs[,-(2:3)]
  fezf2_dep_tras_mTECs$source <- "fezf2_dep_tras_mTECs"
  # other tras in mTECs
  other_tras <- merge(tra_genes, aire_dep_genes, by="entrezgene", all.x = TRUE)
  other_tras <- other_tras[is.na(other_tras$ensembl_gene_id.y),]
  other_tras <- merge(other_tras, fezf2_dep_genes, ny="entrezgene", all.x = TRUE)
  other_tras <- other_tras[is.na(other_tras$ensembl_gene_id),]
  other_tras_mTECs <- merge(other_tras, features_dfs[[i]][[1]], by.x="entrezgene", by.y="geneID")
  other_tras_mTECs <- other_tras_mTECs[,-(2:4)]
  other_tras_mTECs$source <- "other_tras_mTECs"
  # HK in mTECs
  hk_mTECs <- merge(hk_genes, features_dfs[[i]][[1]], by.x="entrezgene", by.y="geneID")
  hk_mTECs <- hk_mTECs[,-2]
  hk_mTECs$source <- "hk_mTECs"
  # mTEC specific in mTECs
  mTEC_spec_mTECs <- merge(mTEC_genes, features_dfs[[i]][[1]], by.x="entrezgene", by.y="geneID")
  mTEC_spec_mTECs$source <- "mTEC_spec_mTECs"
  # thymus specific in mTECs
  thymus_spec_mTECs <- merge(thymus_genes, features_dfs[[i]][[1]], by.x="entrezgene", by.y="geneID")
  thymus_spec_mTECs$source <- "thymus_spec_mTECs"
  # HK in tissues
  hk_tissues <- merge(hk_genes, features_dfs[[i]][[2]], by.x="entrezgene", by.y="geneID")
  hk_tissues <- hk_tissues[,-2]
  hk_tissues$source <- "hk_tissues"
  # tras in tissues 
  tra_tissues <- merge(tra_genes, features_dfs[[i]][[2]], by.x="entrezgene", by.y="geneID")
  tra_tissues <- tra_tissues[,-2]
  tra_tissues$source <- "tra_tissues"
  
  dfs <- list(aire_dep_tras_mTECs, fezf2_dep_tras_mTECs, other_tras_mTECs, hk_mTECs, mTEC_spec_mTECs, thymus_spec_mTECs, hk_tissues, tra_tissues)
  
  plot_dfs <- list()
  c <- 1
  pdf(file = paste("/home/stroemic/hiwi_16/analysis/gene_lists/plots/absolut_positions_percentage_",class,".pdf", sep=""))
  for (i in c(10,11,14,17,15,16,21,22,23)) {
    name <- colnames(other_tras_mTECs)[i]
    print(name)
    ##define df for plotting using percentage function
    plot_df <- data.frame("feature"=sapply(dfs, function(x) percentage(x[!grepl("inter", x$entrezgene),],i)), 
                          "source"=factor(c("aire_dep_tras_mTECs", "fezf2_dep_tras_mTECs", "other_tras_mTECs", "hk_mTECs", "mTEC_spec_mTECs", "thymus_spec_mTECs", "hk_tissues", "tra_tissues"), levels = c("aire_dep_tras_mTECs", "fezf2_dep_tras_mTECs", "other_tras_mTECs", "hk_mTECs", "mTEC_spec_mTECs", "thymus_spec_mTECs", "hk_tissues", "tra_tissues")),
                          "iteration"=name)
    plot_dfs[[c]] <- plot_df
    c=c+1
    print(ggplot(plot_df, aes(x=source, y=feature)) + 
            geom_bar(stat = "identity", fill="#2171b5") +
            labs(title=paste("Feature \"",name,"\" in different datasets", sep = ""), x ="", y =name ) +
            theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
            scale_y_continuous(labels = percent)
    )
  }
  dev.off()
  
  plot_dfs_all <- do.call(rbind, plot_dfs[c(1:6)])
  print(unique(plot_dfs_all$iteration))
  pdf(file = paste("/home/stroemic/hiwi_16/analysis/gene_lists/plots/absolut_positions_proportions_",class,".pdf", sep=""))
  print(ggplot(plot_dfs_all, aes(x=source,y=feature,fill=iteration)) + 
          geom_bar(stat = "identity") +
          labs(title="Percentage of tags belonging to each feature in different datasets", x ="", y ="Percentage") +
          theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
          scale_fill_manual(values = c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')) +
          scale_y_continuous(labels = percent) +
          guides(fill=guide_legend(title=NULL))

  )
  dev.off()
  
  
  
}