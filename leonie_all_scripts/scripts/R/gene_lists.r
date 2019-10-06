library(dplyr)
library(biomaRt)


### read in data frames 
tra_human <- read.table("/home/stroemic/hiwi_16/documentation/gene_lists/tra.human.index2.5x.genes.annotated.tissues.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
tra_mm9_san <- read.table("/home/stroemic/hiwi_16/documentation/gene_lists/sansom_supp.2_tra.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
### aire dep old list
#aire_dep <- read.table("/home/stroemic/hiwi_16/documentation/gene_lists/Aire_dep_hi.csv", header = TRUE, sep = ",")
### aire dep new list sansom et al.
aire_dep_san <- read.table("/home/stroemic/hiwi_16/documentation/gene_lists/sansom_supp.3_table16.csv", header = TRUE, sep = ",")
#all_tissues <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.filtered.csv", header = TRUE, sep = ",")
#all_tissues_wo <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.wo.thymus.filtered.csv", header = TRUE, sep = ",")
#mTECs <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered.csv", header = TRUE, sep = ",")
all_tissues_expression <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/gene_expression_stats_all_tissues.csv", header = TRUE, sep = ",")
housekeeping <- read.table("/home/stroemic/hiwi_16/documentation/gene_lists/housekeeping_genes.tsv", header = TRUE, sep = "\t")
fezf2 <- read.table("/home/stroemic/hiwi_16/documentation/gene_lists/140417mTECmicroarray.csv", header = TRUE, sep = ",", na.strings = "")

#watch out!! other genes in all_tissues_expression than in all_tissues and all_tissues_wo because of differing filtering with rf!!!!

### define mart to map ensembl to entrez IDs
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="grch37.ensembl.org", 
                   dataset="hsapiens_gene_ensembl")
mm_ensembl<- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                       dataset="mmusculus_gene_ensembl")

### define list of TRAs from Sansom
tra_mm9_san <- tra_mm9_san[tra_mm9_san$GNF_GeneAtlas_Specificity == "TRE",]
tras_mm9 <- getBM(attributes =c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'), 
                     filters = 'ensembl_gene_id', 
                     values = tra_mm9_san$Ensembl_ID, 
                     mart = mm_ensembl)
tras_human <- getBM(attributes =c('ensembl_gene_id','entrezgene'), 
                            filters = 'ensembl_gene_id', 
                            values = tras_mm9$hsapiens_homolog_ensembl_gene, 
                            mart = ensembl)
tras_human <- tras_human[!duplicated(tras_human$entrezgene),]
tras_human <- tras_human[!is.na(tras_human$entrezgene),]

# ### define list with TRAs, old list
# tra_entrez <- data.frame("ensembl_gene_id"=tra_human$gene.ids,"entrezgene"=tra_human$entrezID)
# tra_entrez <- tra_entrez[!is.na(tra_entrez$entrezgene),]
# tra_entrez$split <- strsplit(as.character(tra_entrez$entrezgene), "/")
# tra_entrez$entrezgene <- as.character(sapply(tra_entrez$split, function(x) min(as.numeric(x))))
# tra_entrez <- tra_entrez[1:2]
# tra_entrez <- tra_entrez[!(duplicated(tra_entrez$entrezgene)),]

### define list with aire dependent genes new list, sansom et al.
aire_dep_san <- aire_dep_san[(aire_dep_san$fold_change <= 2) &(aire_dep_san$FDR_BH < 0.05),]
mm_aire_dep <- getBM(attributes =c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'), 
                      filters = 'ensembl_gene_id', 
                      values = aire_dep_san$Ensembl, 
                      mart = mm_ensembl)
human_aire_dep_san <- getBM(attributes =c('ensembl_gene_id','entrezgene'), 
                         filters = 'ensembl_gene_id', 
                         values = mm_aire_dep$hsapiens_homolog_ensembl_gene, 
                         mart = ensembl)
human_aire_dep_san <- human_aire_dep_san[!duplicated(human_aire_dep_san$entrezgene),]

# ### define list with aire dependent genes, old list
# aire_dep <- aire_dep[!is.na(aire_dep$human.homologs),]
# human_aire_dep <- getBM(attributes =c('ensembl_gene_id','entrezgene'), 
#                         filters = 'ensembl_gene_id', 
#                         values = aire_dep$human.homologs, 
#                         mart = ensembl)
# human_aire_dep <- human_aire_dep[!(duplicated(human_aire_dep$entrezgene)),]

### define list with Fezf2 dependent genes
fezf2_mm9 <- fezf2[(fezf2$X <= -2)&(fezf2$Student.TTEST < 0.05),]
fezf2_mm9 <- fezf2_mm9[!is.na(fezf2_mm9$Entrez.Gene),]
fezf2_mm9$split <- strsplit(as.character(fezf2_mm9$Entrez.Gene), "///")
fezf2_mm9$entrezgene <- as.character(sapply(fezf2_mm9$split, function(x) min(as.numeric(x))))
fezf2_mm9 <- fezf2_mm9[!duplicated(fezf2_mm9$entrezgene),]
mm_fezf2_dep <- getBM(attributes =c('ensembl_gene_id','entrezgene'), 
                      filters = 'entrezgene', 
                      values = fezf2_mm9$entrezgene, 
                      mart = mm_ensembl)
mm_fezf2_dep_hu <- getBM(attributes =c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'), 
                         filters = 'ensembl_gene_id', 
                         values = mm_fezf2_dep$ensembl_gene_id, 
                         mart = mm_ensembl)
human_fezf2_dep <- getBM(attributes =c('ensembl_gene_id','entrezgene'), 
                          filters = 'ensembl_gene_id', 
                          values = mm_fezf2_dep_hu$hsapiens_homolog_ensembl_gene, 
                          mart = ensembl)
human_fezf2_dep <- human_fezf2_dep[!duplicated(human_fezf2_dep$entrezgene),]

### define list with thymus specific genes
thymus_specific <- all_tissues_expression[(all_tissues_expression$thymus > (10*all_tissues_expression$mean)),c(45,55:60)]
testis_specific <- all_tissues_expression[(all_tissues_expression$testis > (10*all_tissues_expression$mean)),c(43,55:60)]
mTEC_specific <- all_tissues_expression[(all_tissues_expression$all_mTECs > (10*all_tissues_expression$mean)),c(54,55:60)]
thymus_specific_genes <- data.frame("entrezgene"=row.names(thymus_specific))
testis_specific_genes <- data.frame("entrezgene"=row.names(testis_specific))
mTEC_specific_genes <- data.frame("entrezgene"=row.names(mTEC_specific))


###define housekeeping genes 
housekeeping_genes <- getBM(attributes =c('refseq_mrna','entrezgene'), 
                        filters = 'refseq_mrna', 
                        values = housekeeping$refseq, 
                        mart = ensembl)
housekeeping_genes <- housekeeping_genes[!(duplicated(housekeeping_genes$entrezgene)),]
housekeeping_genes <- housekeeping_genes[!is.na(housekeeping_genes$entrezgene),]

##table with tras
write.table(tra_entrez,file ="/home/stroemic/hiwi_16/analysis/gene_lists/human_tra_genes.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(tras_human,file ="/home/stroemic/hiwi_16/analysis/gene_lists/human_tra_genes_sansom.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
##table with aire dep genes 
write.table(human_aire_dep,file ="/home/stroemic/hiwi_16/analysis/gene_lists/human_aire_dep_genes.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(human_aire_dep_san,file ="/home/stroemic/hiwi_16/analysis/gene_lists/human_aire_dep_san_genes.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

##table with aire dep genes 
write.table(human_fezf2_dep,file ="/home/stroemic/hiwi_16/analysis/gene_lists/human_fezf2_dep_genes.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
##table with thymus specific genes 
write.table(thymus_specific_genes,file ="/home/stroemic/hiwi_16/analysis/gene_lists/thymus_specific_genes.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
##table with testis specific genes
write.table(testis_specific_genes,file ="/home/stroemic/hiwi_16/analysis/gene_lists/testis_specific_genes.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
##table with mTEC specific genes
write.table(mTEC_specific_genes,file ="/home/stroemic/hiwi_16/analysis/gene_lists/mTEC_specific_genes.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
##table with housekeeping genes defined on all tissues
write.table(housekeeping_genes,file ="/home/stroemic/hiwi_16/analysis/gene_lists/housekeeping_genes.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")




#########################################################
# all_tissues_genes <- all_tissues[!duplicated(all_tissues[,c("geneID")]),c(2,6)]
# all_tissues_wo_genes <- all_tissues_wo[!duplicated(all_tissues_wo[,c("geneID")]),c(2,6)]
# mTECs_genes <- mTECs[!duplicated(mTECs[,c("geneID")]),c(2,6)]
# 
# thymus_specific <- merge(all_tissues_genes, all_tissues_wo_genes, by=c("geneID"), all.x=TRUE)
# thymus_specific <- thymus_specific[is.na(thymus_specific$BarcodeCount.y),]
# 
# mTEC_specific <- merge(mTECs_genes, all_tissues_wo_genes, by = c("geneID"), all.x = TRUE)
# mTEC_specific <- mTEC_specific[is.na(mTEC_specific$BarcodeCount.y),]
# mTEC_specific_genes <- getBM(attributes =c('ensembl_gene_id','entrezgene'), 
#                              filters = 'entrezgene', 
#                              values = mTEC_specific$geneID, 
#                              mart = ensembl)

# all_tissues_expression$housekeeping <- apply(all_tissues_expression[,-(54:60)], 1, function(x) length(which(x > 100)) >= (0.9*53))
# housekeeping <- all_tissues_expression[all_tissues_expression$housekeeping == TRUE, c(55:61)]
# housekeeping_genes <- data.frame("entrezgene"=row.names(housekeeping))
