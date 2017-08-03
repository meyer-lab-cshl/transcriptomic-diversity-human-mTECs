####Tests
tx_genes_sub <- genes(txdb, filter=list(cds_chrom="chr2"))
transcripts_genes = transcriptsBy(txdb, by="gene")

####Define GRanges
#TSS +-10
TSS10_genes <- promoters(tx_genes, upstream = 10, downstream = 10)
#TSS +-100
TSS100_genes <- promoters(tx_genes, upstream = 100, downstream = 100)
#First Exon 
first_exon = exonsBy(txdb, by='gene')
head(first_exon[[1]], n=1)
#



#### FindOverlaps
df$overlaps=findOverlaps(df_ranges, tx_genes, select = 'first')
df$regions= with(df, names(tx_genes)[overlaps])
lookup <- unique(anno[c('entrezgene', 'ensembl_gene_id')])
df_new <- merge(x = df, y = lookup, by.x = "regions", by.y = "entrezgene", all.x = TRUE, incomparables = NA)
df_new <- df_new[!duplicated(df_new[,c("name")]),]