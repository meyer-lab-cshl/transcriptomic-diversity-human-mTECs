library(ggplot2)
library(raster)
library(reshape2)
library(scales)
library(RColorBrewer)

##read in count matrices
counts_mTECs <- readRDS(file="/home/stroemic/hiwi_16/analysis/gene_lists/main_tss_new_classes_tss_oe_tra_san_aire_san_counts_mTECs_matrix_human_filt.rds")
counts_tissues <- readRDS(file="/home/stroemic/hiwi_16/analysis/gene_lists/main_tss_new_classes_tss_oe_tra_san_aire_san_counts_tissues_matrix_human_filt.rds")


#### further calculations and plotting ####
counts <- cbind(counts_mTECs, counts_tissues)
counts_trunc <- counts[-14,-c(6,12)]
counts_trunc <- t(t(counts_trunc)/rowSums(t(counts_trunc)))

####differing analysis####
###1. Main TSS pattern
#A comparison mTECs / tissues
counts_TSS <- counts[c(1:2),c(1:4,7:10)] ####change here
m_counts_TSS <- melt(counts_TSS)
m_counts_TSS <- transform(m_counts_TSS, feature_class = factor(Var1, rownames(counts_TSS)),
                      gene_class = factor(Var2, colnames(counts_TSS)))

pdf(file ="/home/stroemic/hiwi_16/analysis/gene_lists/plots/Final/main_TSS_patterns_mTECs_vs_tissues_final.pdf")
my.labs <- list("main TSS \u00B1 10","main TSS \u00B1 100")
my.x.labs <- list("Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping","Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping")
print(ggplot(m_counts_TSS, aes(x=gene_class,y=value,fill=feature_class)) + 
        geom_bar(position = "fill", stat = "identity") +
        labs(title="Proportion of tags belonging to each feature in different gene sets", x ="", y ="Proportion") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
        scale_x_discrete(labels = my.x.labs) +
        #scale_fill_brewer(palette="Blues") +
        scale_fill_manual(values = c('#abd9e9','#74add1'), labels = my.labs) +
        guides(fill=guide_legend(title=NULL))
)
dev.off()

#B mTECs only
counts_TSS <- counts[c(1:2),c(1:4)] ####change here
m_counts_TSS <- melt(counts_TSS)
m_counts_TSS <- transform(m_counts_TSS, feature_class = factor(Var1, rownames(counts_TSS)),
                          gene_class = factor(Var2, colnames(counts_TSS)))

pdf(file ="/home/stroemic/hiwi_16/analysis/gene_lists/plots/Final/main_TSS_patterns_mTECs_only_final.pdf")
my.labs <- list("main TSS \u00B1 10","main TSS \u00B1 100")
my.x.labs <- list("Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping")
print(ggplot(m_counts_TSS, aes(x=gene_class,y=value,fill=feature_class)) + 
        geom_bar(position = "fill", stat = "identity") +
        labs(title="Proportion of tags belonging to each feature in different gene sets", x ="", y ="Proportion") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
        scale_x_discrete(labels = my.x.labs) +
        #scale_fill_brewer(palette="Blues") +
        scale_fill_manual(values = c('#abd9e9','#74add1'), labels = my.labs) +
        guides(fill=guide_legend(title=NULL))
      
)
dev.off()




###2. Canonical TSS pattern 
#A comparison mTECs / tissues
perc_TSS <- counts_trunc[c(3:5),c(1:4,6:9)]
m_perc_TSS <- melt(perc_TSS)
m_perc_TSS <- transform(m_perc_TSS, feature_class = factor(Var1, rownames(perc_TSS)),
                      gene_class = factor(Var2, colnames(perc_TSS)))

pdf(file ="/home/stroemic/hiwi_16/analysis/gene_lists/plots/Final/canonical_TSS_patterns_mTECs_vs_tissues_final.pdf")
my.labs <- list("TSS \u00B1 10","TSS \u00B1 100","First exon")
my.x.labs <- list("Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping","Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping")
print(ggplot(m_perc_TSS, aes(x=gene_class,y=value,fill=feature_class)) + 
        geom_bar(stat = "identity") +
        labs(title="Percentage of tags belonging to each feature in different gene sets", x ="", y ="Percentage") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
        scale_x_discrete(labels = my.x.labs) +
        #scale_fill_brewer(palette="Blues") +
        scale_fill_manual(values = c('#abd9e9','#74add1','#4575b4'), labels = my.labs) +
        guides(fill=guide_legend(title=NULL)) +
		scale_y_continuous(labels = percent)
)
dev.off()
#B mTECs only
perc_TSS <- counts_trunc[c(3:5),c(1:4)] ####change here
m_perc_TSS <- melt(perc_TSS)
m_perc_TSS <- transform(m_perc_TSS, feature_class = factor(Var1, rownames(perc_TSS)),
                          gene_class = factor(Var2, colnames(perc_TSS)))

pdf(file ="/home/stroemic/hiwi_16/analysis/gene_lists/plots/Final/canonical_TSS_patterns_mTECs_only_final.pdf")
my.labs <- list("TSS \u00B1 10","TSS \u00B1 100","First exon")
my.x.labs <- list("Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping")
print(ggplot(m_perc_TSS, aes(x=gene_class,y=value,fill=feature_class)) + 
        geom_bar(stat = "identity") +
        labs(title="Percentage of tags belonging to each feature in different gene sets", x ="", y ="Percentage") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
        scale_x_discrete(labels = my.x.labs) +
        #scale_fill_brewer(palette="Blues") +
        scale_fill_manual(values = c('#abd9e9','#74add1','#4575b4'), labels = my.labs) +
        guides(fill=guide_legend(title=NULL)) +
		scale_y_continuous(labels = percent)
      
)
dev.off()


###2. Intronic/Antisense high in Aire regulated TRAs
#A comparison mTECs / tissues
perc_intron_anti <- counts_trunc[c(11:12),c(1:4,6:9)] ####change here
m_perc_intron_anti <- melt(perc_intron_anti)
m_perc_intron_anti <- transform(m_perc_intron_anti, feature_class = factor(Var1, rownames(perc_intron_anti)),
                          gene_class = factor(Var2, colnames(perc_intron_anti)))

pdf(file ="/home/stroemic/hiwi_16/analysis/gene_lists/plots/Final/intronic_antisense_mTECs_tissues_with_main_final.pdf")
my.labs <- list("Intronic","Antisense")
my.x.labs <- list("Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping","Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping")
print(ggplot(m_perc_intron_anti, aes(x=gene_class,y=value,fill=feature_class)) + 
        geom_bar(stat = "identity") +
        labs(title="Proportion of tags belonging to each feature in different gene sets", x ="", y ="Proportion") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
        #scale_fill_brewer(palette="Blues") +
        scale_fill_manual(values = c('#fdae61','#d73027'), labels = my.labs) +
        scale_x_discrete(labels = my.x.labs) +
        guides(fill=guide_legend(title=NULL)) +
        scale_y_continuous(labels = percent)
      
)
dev.off()

#B comparison mTECs / tissues with other protein coding genes 
perc_intron_anti <- counts_trunc[c(11:12),c(1:5,6:10)] ####change here
m_perc_intron_anti <- melt(perc_intron_anti)
m_perc_intron_anti <- transform(m_perc_intron_anti, feature_class = factor(Var1, rownames(perc_intron_anti)),
                          gene_class = factor(Var2, colnames(perc_intron_anti)))

pdf(file ="/home/stroemic/hiwi_16/analysis/gene_lists/plots/Final/intronic_antisense_mTECs_tissues_with_main_final_protein_coding_genes.pdf")
my.labs <- list("Intronic","Antisense")
my.x.labs <- list("Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping","Other protein\ncoding genes","Aire TRAs","Fezf2 TRAs","Other TRAs","Housekeeping","Other protein\ncoding genes")
print(ggplot(m_perc_intron_anti, aes(x=gene_class,y=value,fill=feature_class)) + 
        geom_bar(stat = "identity") +
        labs(title="Proportion of tags belonging to each feature in different gene sets", x ="", y ="Proportion") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
        #scale_fill_brewer(palette="Blues") +
        scale_fill_manual(values = c('#fdae61','#d73027'), labels = my.labs) +
        scale_x_discrete(labels = my.x.labs) +
        guides(fill=guide_legend(title=NULL)) +
        scale_y_continuous(labels = percent)
      
)
dev.off()


###3. TSS at other exons 
#A comparison mTECs / tissues as whole with newly calculated percentages --> real data 
df_tss_oe <- data.frame("feature"=rownames(counts), "tras_mTECs"=(counts[,1]+counts[,2]+counts[,3]),
                        "hk_mTECs"=counts[,4],
                        "tras_tissues"=(counts[,7]+counts[,8]+counts[,9]),
                        "hk_tissues"=counts[,10])

tss_oe <- data.matrix(df_tss_oe[,2:5])
perc_tss_oe <- t(t(tss_oe)/rowSums(t(tss_oe)))
perc_tss_oe <- perc_tss_oe[6:7,]
m_perc_tss_oe <- melt(perc_tss_oe)
m_perc_tss_oe<- transform(m_perc_tss_oe, feature_class = factor(Var1, rownames(perc_tss_oe)),
                          gene_class = factor(Var2, colnames(perc_tss_oe)))

pdf(file ="/home/stroemic/hiwi_16/analysis/gene_lists/plots/Final/tss_oe_mTECs_tissues_with_main_final.pdf")
my.labs <- list("Downstream exon start \u00B1 10","Downstream exon start \u00B1 100")
my.x.labs <- list("TRAs","Housekeeping","TRAs","Housekeeping")
print(ggplot(m_perc_tss_oe, aes(x=gene_class,y=value,fill=feature_class)) + 
        geom_bar(stat = "identity") +
        labs(title="Proportion of tags belonging to each feature in different gene sets", x ="", y ="Proportion") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
        #scale_fill_brewer(palette="Blues") +
        scale_fill_manual(values = c('#a6d96a','#1a9641'), labels = my.labs) +
        scale_x_discrete(labels = my.x.labs) +
        guides(fill=guide_legend(title=NULL)) +
        scale_y_continuous(labels = percent)
      
)
dev.off()


#B comparison mTECs / tissues as whole with newly calculated percentages --> real data with protein coding
df_tss_oe <- data.frame("feature"=rownames(counts), "tras_mTECs"=(counts[,1]+counts[,2]+counts[,3]),
                        "hk_mTECs"=counts[,4], 
                        "others_mTECs"=counts[,5],
                        "tras_tissues"=(counts[,7]+counts[,8]+counts[,9]),
                        "hk_tissues"=counts[,10],
                        "others_tissues"=counts[,11])

tss_oe <- data.matrix(df_tss_oe[,2:7])
perc_tss_oe <- t(t(tss_oe)/rowSums(t(tss_oe)))
perc_tss_oe <- perc_tss_oe[6:7,]
m_perc_tss_oe <- melt(perc_tss_oe)
m_perc_tss_oe<- transform(m_perc_tss_oe, feature_class = factor(Var1, rownames(perc_tss_oe)),
                          gene_class = factor(Var2, colnames(perc_tss_oe)))

pdf(file ="/home/stroemic/hiwi_16/analysis/gene_lists/plots/Final/tss_oe_mTECs_tissues_with_main_final_protein_coding_genes.pdf")
my.labs <- list("Downstream exon start \u00B1 10","Downstream exon start \u00B1 100")
my.x.labs <- list("TRAs","Housekeeping","Other protein\ncoding genes","TRAs","Housekeeping","Other protein\ncoding genes")
print(ggplot(m_perc_tss_oe, aes(x=gene_class,y=value,fill=feature_class)) + 
         geom_bar(stat = "identity") +
        labs(title="Proportion of tags belonging to each feature in different gene sets", x ="", y ="Proportion") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=10)) +
        #scale_fill_brewer(palette="Blues") +
        scale_fill_manual(values = c('#a6d96a','#1a9641'), labels = my.labs) +
        scale_x_discrete(labels = my.x.labs) +
        guides(fill=guide_legend(title=NULL)) +
        scale_y_continuous(labels = percent)
      
)
dev.off()
