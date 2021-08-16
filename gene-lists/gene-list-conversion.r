## libraries ####
library(tidyverse)
library(biomaRt)

indir <- "~/data/tss"

## read in gene lists from literature ####
## aire dependent genes from Sansom et a ().Supp Tab 3, sheet 16 ####
aire<- read_csv(file.path(indir, "gene_lists/input/sansom_supp.3_table16.csv"))

## housekeeping genes from Eisenberg & Levanon 2013
housekeeping <- read_tsv(file.path(indir, "gene_lists/input/housekeeping_genes.tsv"))

## fezf2 genes from Takaba 2015
fezf2 <- read_csv(file.path(indir, "gene_lists/input/140417mTECmicroarray.csv"))

## define marts to map mouse to human and ensembl to entrez IDs ####
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl")
mm_ensembl<- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="mmusculus_gene_ensembl")

## define list of human aire dependent genes based on sansom et al. ####
## Fold change column defined by authors as aire_ko/aire_wt,
## ie small fold changes indicate upregulation through aire
aire <- aire[(aire$fold_change <= 2) & (aire$FDR_BH < 0.05),]
write_csv(aire,
          file.path(indir, "gene_lists/input/mm_aire-dependent_sansom_supp3_table16.csv"))
mm_aire <- getBM(attributes =c('ensembl_gene_id',
                               'hsapiens_homolog_ensembl_gene'),
                      filters = 'ensembl_gene_id',
                      values = aire$Ensembl,
                      mart = mm_ensembl)
human_aire <- getBM(attributes =c('ensembl_gene_id','entrezgene_id'),
                         filters = 'ensembl_gene_id',
                         values = mm_aire$hsapiens_homolog_ensembl_gene,
                         mart = ensembl)
human_aire <- human_aire[!duplicated(human_aire$entrezgene),]

## define list with Fezf2 dependent genes based on Takaba 2015 ####
## Fold change column defined by authors as -1* WT/KO,
## ie large negative fold changes indicate upregulation through fezf2
fezf2_mm9 <- fezf2 %>%
    filter(FoldChange <= -2, StudentTTEST < 0.05) %>%
    drop_na(EntrezGene) %>%
    separate_rows(EntrezGene, sep="///", convert = TRUE)
fezf2_mm9 <- fezf2_mm9[!duplicated(fezf2_mm9$EntrezGene),]

mm_fezf2 <- getBM(attributes =c('ensembl_gene_id','entrezgene_id'),
                      filters = 'entrezgene_id',
                      values = fezf2_mm9$EntrezGene,
                      mart = mm_ensembl)
fezf2_mm9 <- fezf2_mm9 %>%
    left_join(mm_fezf2, by=c("EntrezGene"="entrezgene_id")) %>%
    dplyr::select(ensembl_gene_id, EntrezGene, everything()) %>%
    rename(Ensembl="ensembl_gene_id")
write_csv(fezf2_mm9,
          file.path(indir, "gene_lists/input/mm_fezf2-dependent_takaba_140417mTECmicroarray.csv"))


mm_fezf2_hu <- getBM(attributes =c('ensembl_gene_id',
                                       'hsapiens_homolog_ensembl_gene'),
                         filters = 'ensembl_gene_id',
                         values = mm_fezf2$ensembl_gene_id,
                         mart = mm_ensembl)
human_fezf2 <- getBM(attributes =c('ensembl_gene_id','entrezgene_id'),
                          filters = 'ensembl_gene_id',
                          values = mm_fezf2_hu$hsapiens_homolog_ensembl_gene,
                          mart = ensembl)
human_fezf2 <- human_fezf2[!duplicated(human_fezf2$entrezgene_id),]

## define housekeeping genes ####
housekeeping_genes <- getBM(attributes =c('refseq_mrna','entrezgene_id'),
                        filters = 'refseq_mrna',
                        values = housekeeping$refseq,
                        mart = ensembl)
housekeeping_genes <- housekeeping_genes[!(duplicated(housekeeping_genes$entrezgene)),]

## write tables with gene lists ####
write_csv(human_aire,
          file.path(indir, "gene_lists/human_aire_dep_genes_san.csv"))
write_csv(human_fezf2,
          file.path(indir, "gene_lists/human_fezf2_dep_genes.csv"))
write.table(housekeeping_genes,
            file.path(indir, "gene_lists/housekeeping_genes.csv"))

