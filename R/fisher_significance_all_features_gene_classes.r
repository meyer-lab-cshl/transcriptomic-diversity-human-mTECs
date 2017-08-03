
counts_mTECs <- readRDS(file="/home/stroemic/hiwi_16/analysis/gene_lists/tss_oe_tra_san_aire_san_counts_mTECs_matrix_human_filt.rds")
counts_tissues <- readRDS(file="/home/stroemic/hiwi_16/analysis/gene_lists/tss_oe_tra_san_aire_san_counts_tissues_matrix_human_filt.rds")
counts <- cbind(counts_mTECs, counts_tissues)


counts_2 <- rbind(counts[1,], t(rowSums(t(counts[2:11,]))))
fisher.test(counts_2[,c(1:3)])