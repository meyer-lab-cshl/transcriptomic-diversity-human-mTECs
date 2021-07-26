## libraries ####
library(tidyverse)
library(biomaRt)

indir <- "~/data/tss/gene_lists/public_data"
tradir <- "~/data/tra"
outdir <- "~/data/tss/gene_lists/gene_lists_human"

## read in gene lists from literature ####
## aire dependent genes from Sansom et a ().Supp Tab 3, sheet 16 ####
aire<- read_csv(file.path(indir, "sansom_supp.3_table16.csv"))

## housekeeping genes from Eisenberg & Levanon 2013
housekeeping <- read_tsv(file.path(indir, "housekeeping_genes.tsv"))

## fezf2 genes from Takaba 2015
fezf2 <- read_csv(file.path(indir, "140417mTECmicroarray.csv"),
                  col_names = c("ProveSetID", "WTaverage", "KOaverage",
                                "FoldChange", "StudentTTEST", "GeneSymbol",
                                "EntrezGene"),
                  skip=1
                  )

## TRA genes and transcripts from tra/tra.R ####
tra_genes <- read_csv(file.path(tradir,
                                "gtex_genes_binarized_expression.csv"))
tra_transcripts <- read_csv(file.path(tradir,
                                      "gtex_transcripts_binarized_expression.csv"))

## define marts to map mouse to human and ensembl to entrez IDs ####
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl")
mm_ensembl<- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="mmusculus_gene_ensembl")

## define list of human aire dependent genes based on sansom et al. ####
## Fold change column defined by authors as aire_ko/aire_wt,
## ie small fold changes indicate upregulation through aire
aire <- aire[(aire$fold_change <= 2) & (aire$FDR_BH < 0.05),]
mm_aire <- getBM(attributes =c('ensembl_gene_id',
                               'hsapiens_homolog_ensembl_gene'),
                      filters = 'ensembl_gene_id',
                      values = aire$Ensembl,
                      mart = mm_ensembl)
human_aire <- getBM(attributes = c('refseq_mrna', 'entrezgene_id',
                                  'ensembl_gene_id', 'external_gene_name'),
                         filters = 'ensembl_gene_id',
                         values = mm_aire$hsapiens_homolog_ensembl_gene,
                         mart = ensembl)
human_aire <- human_aire[!duplicated(human_aire$ensembl_gene_id),]

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
mm_fezf2_hu <- getBM(attributes =c('ensembl_gene_id',
                                       'hsapiens_homolog_ensembl_gene'),
                         filters = 'ensembl_gene_id',
                         values = mm_fezf2$ensembl_gene_id,
                         mart = mm_ensembl)
human_fezf2 <- getBM(attributes =c('refseq_mrna', 'entrezgene_id',
                                   'ensembl_gene_id', 'external_gene_name'),
                          filters = 'ensembl_gene_id',
                          values = mm_fezf2_hu$hsapiens_homolog_ensembl_gene,
                          mart = ensembl)
human_fezf2 <- human_fezf2[!duplicated(human_fezf2$ensembl_gene_id),]

## define housekeeping genes ####
housekeeping_genes <- getBM(attributes =c('refseq_mrna','entrezgene_id',
                                          'ensembl_gene_id',
                                          'external_gene_name'),
                        filters = 'refseq_mrna',
                        values = housekeeping$refseq,
                        mart = ensembl)
housekeeping_genes <- housekeeping_genes[!(duplicated(housekeeping_genes$ensembl_gene_id)),]

## define TRA genes and transcripts ####
tra_genes_annotated <- getBM(attributes =c('refseq_mrna','entrezgene_id',
                                          'ensembl_gene_id',
                                          'external_gene_name'),
                            filters = 'ensembl_gene_id',
                            values = tra_genes$id,
                            mart = ensembl)
tra_genes_annotated <- tra_genes_annotated[!(duplicated(tra_genes_annotated$ensembl_gene_id)),]

tra_genes_annotated <- tra_genes_annotated %>%
    left_join(tra_genes, by=c("ensembl_gene_id"="id"))

tra_transcripts_annotated <- getBM(attributes =c('refseq_mrna','entrezgene_id',
                                 'ensembl_gene_id',
                                 'external_gene_name', 'ensembl_transcript_id'),
                   filters = 'ensembl_transcript_id',
                   values = tra_transcripts$id,
                   mart = ensembl)
tra_transcripts_annotated <- tra_transcripts_annotated[!(duplicated(tra_transcripts_annotated$ensembl_transcript_id)),]

tra_transcripts_annotated <- tra_transcripts_annotated %>%
    left_join(tra_transcripts, by=c("ensembl_transcript_id"="id"))

## Remove any ambiguous genes ####

common_genes <- c(tra_genes_annotated$ensembl_gene_id,
                  housekeeping_genes$ensembl_gene_id)
ambiguous_genes <- table(common_genes)[table(common_genes) == 2]

tra_genes_annotated <- tra_genes_annotated %>%
    dplyr::filter(!ensembl_gene_id %in% names(ambiguous_genes))

tra_transcripts_annotated <- tra_transcripts_annotated %>%
    dplyr::filter(!ensembl_gene_id %in% names(ambiguous_genes))

housekeeping_genes <- housekeeping_genes %>%
    dplyr::filter(!ensembl_gene_id %in% names(ambiguous_genes))


## write tables with gene lists ####
write_csv(human_aire,
          file.path(outdir, "human_aire_dep_genes_san.csv"))
write_csv(human_fezf2,
          file.path(outdir, "human_fezf2_dep_genes.csv"))
write_csv(housekeeping_genes,
            file.path(outdir, "housekeeping_genes.csv"))
write_csv(tra_genes_annotated,
          file.path(outdir, "tra_genes.csv"))
write_csv(tra_transcripts_annotated,
          file.path(outdir, "tra_transcripts.csv"))

