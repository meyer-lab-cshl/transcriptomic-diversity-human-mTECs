library(reshape)
library(LSD)
library(dplyr)
library(ggplot2)
library(scales)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)

#### functions ####
preprocess_df <- function(x, keep) {
  x <- x[x$geneID %in% names(cds),]
  x <- x[x$geneID %in% keep,]
  x <- droplevels(x)
  x <- x[with(x, order(chrom, geneID, start)),]
  filt <- x %>% group_by(geneID) %>% mutate(the_rank  = rank(-BarcodeCount, ties.method = "first")) %>% filter(the_rank == 1)
  
  filt$geneID <- as.character(filt$geneID)
  return(filt)
}

##### BioMart object
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="grch37.ensembl.org", 
                   dataset="hsapiens_gene_ensembl")


##mtec all positions
mTEC_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered.csv", header=TRUE, sep=",")
##tissue (wo thymus) all positions
tissues_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.wo.thymus.filtered.csv", header = TRUE, sep = ",")
##get common genes for calculation
keep <- intersect(mTEC_pos_fil$geneID, tissues_pos_fil$geneID)

#### read in dataframe of genes with shifted TSS pattern
df_genes <- read.table("/home/stroemic/hiwi_16/analysis/shifted_TSS/sorted_protein_coding_genes_shifted_tss.csv", header = TRUE, sep=",")


#### read in txdb object and get CDS GRanges

##### txdb gene object
genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)

###mart object
### keep only protein coding genes 
mart_gene_type <- getBM(attributes = c('entrezgene', 'gene_biotype'), 
                        filters = 'entrezgene',
                        values = names(tx_genes), 
                        mart=ensembl)
mart_gene_type <- mart_gene_type[grepl("protein_coding", mart_gene_type$gene_biotype, ignore.case = TRUE),]

cds <- cdsBy(txdb, by="gene")
cds <- cds[names(cds) %in% mart_gene_type$entrezgene]
cds <- cds[names(cds) %in% keep]
oc <- order(as.numeric(names(cds)))
cds <- cds[oc]



#### analysis ####
mTEC_pos_fil <- preprocess_df(mTEC_pos_fil, keep)
tissues_pos_fil <- preprocess_df(tissues_pos_fil, keep)

list_genes_names <- names(cds)
list_peptides <- vector("list", length(list_genes_names))
names(list_peptides) <- list_genes_names


start.time <- Sys.time()
for (i in names(list_peptides)) {
  print(i)
  
  cds_ranges <- cds[[i]]
  cds_ranges <- reduce(cds_ranges)
  main_peak <- mTEC_pos_fil[mTEC_pos_fil$geneID == i,]$start
  
  if(unique(as.character(strand(cds_ranges))) == "+") { # positive strand
    if(main_peak <= start(cds_ranges)[1]) {
      peptide = "no truncation"
      #print(peptide)
    } else if(main_peak >= end(cds_ranges)[length(cds_ranges)]) {
      peptide = "fully truncated"
      #print(peptide)
    } else {
      cds_ranges <- restrict(cds_ranges, start = start(cds_ranges)[1], end = main_peak)
      sequences <- getSeq(genome, cds_ranges)
      sequence <- do.call(paste0, sequences)
      peptide <- translate(DNAString(sequence))
      #print(peptide)
    }
    
  } else { #negative strand
    if(main_peak >= end(cds_ranges)[length(cds_ranges)]) {
      peptide = "no truncation"
      #print(peptide)
    } else if(main_peak <= start(cds_ranges)[1]) {
      ppeptide = "fully truncated"
      #print(peptide)
    } else {
      cds_ranges <- restrict(cds_ranges, start = main_peak, end = end(cds_ranges)[length(cds_ranges)])
      sequences <- getSeq(genome, cds_ranges)
      sequence <- do.call(paste0, rev(sequences))
      peptide <- translate(DNAString(sequence))
      #print(peptide)
    }
  }
  list_peptides[[i]] <- as.character(peptide) 
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

saveRDS(list_peptides, file="/home/stroemic/hiwi_16/analysis/shifted_TSS/list_truncated_peptides_protein_coding_genes.rds")


### make peptide dataframe to merge with gene list
peptides <- data.frame("truncated_peptide"=unlist(list_peptides), "geneID" = names(list_peptides))
df_genes <- merge(df_genes, peptides, by="geneID")

df_genes$truncation_length <- 0
df_genes[!grepl("trunca",df_genes$truncated_peptide),]$truncation_length <- nchar(as.character(df_genes[!grepl("trunca",df_genes$truncated_peptide),]$truncated_peptide))

top_genes <- df_genes[with(df_genes, order(p.value, -ratio, -truncation_length)),]
top_tras <- top_genes[top_genes$tra == "tra",]
write.table(top_genes, file="/home/stroemic/hiwi_16/analysis/shifted_TSS/sorted_protein_coding_genes_shifted_tss_plus_truncation.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "", sep = ",")
write.table(top_tras, file="/home/stroemic/hiwi_16/analysis/shifted_TSS/sorted_tras_protein_coding_genes_shifted_tss_plus_truncation.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "", sep = ",")



pdf(file="/home/stroemic/hiwi_16/analysis/shifted_TSS/plots/all_genes_p.value_vs_truncation.pdf")
gene_ids = c("MUC6","MLANA")
ggplot(df_genes, aes(x=log10(truncation_length), y=-log10(p.value))) +
  geom_point(size=0.5) +
  geom_point(data=df_genes[df_genes$external_gene_name %in% gene_ids,], colour="blue", size=2) +
  #geom_vline(xintercept = 0, linetype=2, colour= "red") +
  #scale_x_continuous(limits =c(-4,6.5), oob = squish) +
  #scale_y_continuous(limits=c(0,340),oob=squish) +
  #geom_text(data=subset(df_genes, log(ratio) > 4.8 | log(ratio) < -2.5),aes(label=row.names(subset(df_genes, log(ratio) > 4.8 | log(ratio) < -2.5))), color="red",hjust=0.5,vjust=1) +
  labs(title="Genewise significance of shifted expression pattern between tissues and mTECs \nover truncated gene expression in mTECs", x ="log10(length of protein truncation to main TSS in mTECs)", y ="-log10(p.value)") +
  #annotate("text", x=2.5, y=370, label=paste(nrow(df_genes[df_genes$p.value < 0.05 & df_genes$ratio > 1, ])," sig. genes",sep=""), fontface =2) + 
  #annotate("text", x=-2.5, y=370, label =paste(nrow(df_genes[df_genes$p.value < 0.05 & df_genes$ratio < 1, ])," sig. genes",sep=""), fontface =2) +
  theme_bw()

dev.off()
