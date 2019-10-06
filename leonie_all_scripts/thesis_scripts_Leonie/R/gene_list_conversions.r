library(biomaRt)

##table with thymus specific genes 
thymus_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/thymus_specific_genes.csv", header=TRUE, sep=",")
thymus_genes <- data.frame("entrezgene" =thymus_genes[!grepl("inter", thymus_genes$entrezgene),])
##table with mTEC specific genes
mTEC_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/mTEC_specific_genes.csv", header=TRUE, sep=",")
mTEC_genes <- data.frame("entrezgene" = mTEC_genes[!grepl("inter", mTEC_genes$entrezgene),])


ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="grch37.ensembl.org", 
                   dataset="hsapiens_gene_ensembl")

genes_thymus <- getBM(attributes = c('entrezgene', 'ensembl_gene_id','external_gene_name'), 
               filters = 'entrezgene',
               values = thymus_genes$entrezgene, 
               mart=ensembl)
genes_thymus <- genes_thymus[!duplicated(genes_thymus$entrezgene),]

thymus_tab <- merge(thymus_genes, genes_thymus, by="entrezgene", all.x=TRUE)

genes_mTEC <- getBM(attributes = c('entrezgene', 'ensembl_gene_id','external_gene_name'), 
                      filters = 'entrezgene',
                      values = mTEC_genes$entrezgene, 
                      mart=ensembl)
genes_mTEC <- genes_mTEC[!duplicated(genes_mTEC$entrezgene),]

mTEC_tab <- merge(mTEC_genes, genes_mTEC, by="entrezgene", all.x=TRUE)

write.table(thymus_tab,file ="/home/stroemic/hiwi_16/analysis/gene_lists/thymus_specific_genes_with_names.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(mTEC_tab,file ="/home/stroemic/hiwi_16/analysis/gene_lists/mTEC_specific_genes_with_names.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
