#!/usr/bin/env Rscript

library(sleuth)
library(biomaRt)

base_dir <- '/grid/meyer/home/jacarter/TSS_RNAseq'
kal <- 'Kallisto/Quant'

sample_id <- dir(file.path(base_dir,kal))[2:7]
kal_dirs <- file.path(base_dir,kal,sample_id)

s2c <- read.table(file.path(base_dir,"Code","TSS_sleuth.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~subject+condition, 'full')
so <- sleuth_fit(so, ~subject, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

so <- sleuth_wt(so, which_beta='conditionlow',which_model='full')
sleuth_table <- sleuth_results(so, 'conditionlow', 'wt')
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

write.csv(sleuth_significant, paste(base_dir,'/Sleuth_results_sig.txt',sep=""))
write.csv(sleuth_table, paste(base_dir,'/Sleuth_results_all.txt',sep=""))
