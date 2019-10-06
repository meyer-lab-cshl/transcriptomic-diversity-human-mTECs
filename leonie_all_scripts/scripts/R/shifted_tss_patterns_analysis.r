library(reshape)
library(LSD)
library(ggplot2)
library(scales)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)

##### BioMart object
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="grch37.ensembl.org", 
                   dataset="hsapiens_gene_ensembl")

###read in tra genes 
tra_entrez <- read.table(file ="/home/stroemic/hiwi_16/analysis/gene_lists/human_tra_genes_sansom.csv", sep=",", header = TRUE)
##table with aire dep genes 
aire_dep_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_aire_dep_san_genes.csv", header=TRUE, sep=",")
##table with fezf2 dep genes 
fezf2_dep_genes <- read.table("/home/stroemic/hiwi_16/analysis/gene_lists/human_fezf2_dep_genes.csv", header=TRUE, sep=",")

aire_dep_tras <- merge(tra_entrez, aire_dep_genes, by="entrezgene")
fezf2_dep_tras <- merge(tra_entrez, fezf2_dep_genes, by="entrezgene")


#### read in List of gene matrices
list_genes <- readRDS(file="/home/stroemic/hiwi_16/analysis/shifted_TSS/list_genes_matrices_protein_coding.rds")


#### analysis ####
##calculate p values
list_fisher <- lapply(list_genes, function(x) fisher.test(round(x[,1:2])))

#A) Boxplots of ratios in mTECs and tissues 
##calculate ratios
ratio_tissues<- sapply(list_genes, function(x) x[1,2] / (x[1,1] + x[1,2]))
ratio_mTEC  <- sapply(list_genes, function(x) x[2,2] / (x[2,1] + x[2,2]))
##plot distribution in box_plots
df_ratios <- data.frame("Tissues"=ratio_tissues, "mTECs"=ratio_mTEC)
m_df_ratios <- melt(df_ratios)
pdf(file="/home/stroemic/hiwi_16/analysis/shifted_TSS/plots/boxplot_ratios_tissues_mTECs.pdf")
p1 <- ggplot(m_df_ratios, aes(x=variable, y=value, fill=variable))
ylim1 = boxplot.stats(m_df_ratios$value)$stats[c(1, 5)]
p1 + geom_boxplot(outlier.size=NA) +  coord_cartesian(ylim = ylim1*1.05) + 
    theme_bw() +
    labs(title="Boxplots of ratio between CDS/5'UTR expression \nin tissues and mTECs", x ="", y ="Ratio CDS/(CDS+5'UTR) expression") +
    scale_fill_manual(values = c('#74add1','#4575b4')) +
    guides(fill=FALSE)
dev.off()

#B) Built data frame to filter genes 
#1. Significane
filt_fisher <- sapply(list_fisher, function(x) x$p.value)
filt_fisher <- p.adjust(filt_fisher, method = "BH") 
#2. Shift into CDS
##calculate ratios between tissues and mtecs cds vs. upstream each
filt_ratios <- sapply(list_genes, function(x) (x[2,2]/(x[2,1]+x[2,2]))/(x[1,2] / (x[1,1] + x[1,2])) )
#3. expression levels mTECs, tisssues
tissues_upstream <- sapply(list_genes, function(x) x[1,1])
tissues_cds <- sapply(list_genes, function(x) x[1,2])
mTECs_upstream <- sapply(list_genes, function(x) x[2,1])
mTECs_cds <- sapply(list_genes, function(x) x[2,2])
##geneIDs and values to data.frame 
df <- data.frame("geneID"=as.character(names(filt_ratios)), 
                 "ratio"= filt_ratios, "p.value"=filt_fisher,
                 "tissues_upstream"= tissues_upstream,"tissues_cds"=tissues_cds,
                 "mTECs_upstream"= mTECs_upstream, "mTECs_cds"=mTECs_cds)
##add TRAs
df$tra <- match(df$geneID, tra_entrez$entrezgene)
df[!is.na(df$tra),]$tra<- 'tra'
df[is.na(df$tra),]$tra <- 'non_tra'

df$aire <- match(df$geneID, aire_dep_tras$entrezgene )
df[!is.na(df$aire),]$aire<- 'yes'
df[is.na(df$aire),]$aire <- 'no'

df$fezf2 <- match(df$geneID, fezf2_dep_tras$entrezgene )
df[!is.na(df$fezf2),]$fezf2<- 'yes'
df[is.na(df$fezf2),]$fezf2 <- 'no'

##add gene names
genes <- getBM(attributes = c('entrezgene', 'ensembl_gene_id','external_gene_name'), 
                        filters = 'entrezgene',
                        values = df$geneID, 
                        mart=ensembl)
genes <- genes[!duplicated(genes$entrezgene),]

df_export <- merge(df, genes, by.x="geneID", by.y="entrezgene")

top_genes <- df_export[with(df_export, order(p.value, -ratio)),]
top_tras <- top_genes[top_genes$tra == "tra",]
top_aire_tras <- top_tras[top_tras$aire == "yes",]

write.table(df_export[,c(1,12,11,2:7,8:10)], file="/home/stroemic/hiwi_16/analysis/shifted_TSS/protein_coding_genes_shifted_tss.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "", sep = ",")
write.table(top_genes[,c(1,12,11,2:7,8:10)], file="/home/stroemic/hiwi_16/analysis/shifted_TSS/sorted_protein_coding_genes_shifted_tss.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "", sep = ",")
write.table(top_tras[,c(1,12,11,2:7,8:10)], file="/home/stroemic/hiwi_16/analysis/shifted_TSS/sorted_tras_protein_coding_genes_shifted_tss.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "", sep = ",")
write.table(top_aire_tras[,c(1,12,11,2:7,8:10)], file="/home/stroemic/hiwi_16/analysis/shifted_TSS/sorted_aire_tras_protein_coding_genes_shifted_tss.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "", sep = ",")


####C) plot p.value over ratio #####
# df <- df_export[!is.infinite(df_export$ratio),]
# df <- df[!is.nan(df$ratio),]
# df <- df[df$p.value != 0,]
# df <- df[df$ratio != 0,]


pdf(file="/home/stroemic/hiwi_16/analysis/shifted_TSS/plots/all_genes_p.value_vs_ratio.pdf")
gene_ids = c("MUC6","MLANA","INS")
ggplot(df, aes(x=log(ratio), y=-log10(p.value))) +
        geom_point(size=0.5) +
        geom_point(data=df[df$external_gene_name %in% gene_ids,], colour="blue", size=2) +
        geom_vline(xintercept = 0, linetype=2, colour= "red") +
        #scale_x_continuous(limits =c(-4,6.5), oob = squish) +
        #scale_y_continuous(limits=c(0,340),oob=squish) +
        #geom_text(data=subset(df, log(ratio) > 4.8 | log(ratio) < -2.5),aes(label=row.names(subset(df, log(ratio) > 4.8 | log(ratio) < -2.5))), color="red",hjust=0.5,vjust=1) +
        labs(title="Genewise significance of shifted expression pattern between tissues and mTECs", x ="log(ratio of CDS/5'UTR expression between tissues and mTECs)", y ="-log10(p.value)") +
        annotate("text", x=2.5, y=370, label=paste(nrow(df[df$p.value < 0.05 & df$ratio > 1, ])," sig. genes",sep=""), fontface =2) + 
        annotate("text", x=-2.5, y=370, label =paste(nrow(df[df$p.value < 0.05 & df$ratio < 1, ])," sig. genes",sep=""), fontface =2) +
        theme_bw()

heatscatter(df$ratio, (df$p.value), log = "xy", ylim = rev(range(df$p.value)), 
            main = "Genewise significance of shifted expressio patterns between tissues and mTECs",
            xlab = "log(ratio of CDS/5'UTR expression between tissue and mTECs)",
            ylab = "log(p.value)")
dev.off()

