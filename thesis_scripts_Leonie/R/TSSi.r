library(TSSi)
library(dplyr)
library(plyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(parallel)

setwd("hiwi_16/data/external/Fantom5/CTSS_BED/tissues/")

#####Load annotations for gene conversions
anno <- read.csv("/home/stroemic/hiwi_16/documentation/annotations.csv", header = TRUE, sep = "," )

#####Load Dataset & transform
#df=read.table("/home/stroemic/hiwi_16/data/external/Fantom5/CTSS_BED/tissues/head200_thymus.txt")
fantoms <- list.files(pattern="*.bed.gz")
fantoms_dfs = list()
for(i in fantoms){
  dat <- read.table(i)
  fantoms_dfs[[i]] <- dat
  print(i)
}
df <- bind_rows(fantoms_dfs, .id = NULL)
#rename columns
df<- rename(df,c("V1"="chromosome","V2"="start","V3"="end","V4"="name","V5"="count","V6"="strand"))
df<- subset(df, select = -c(end, name))
df<- aggregate(df$count, by=list(df$chromosome, df$start, df$strand), FUN=sum)
df<- rename(df, c("Group.1"="chromosome", "Group.2"="start", "Group.3"="strand", "x"="count"))

df_plus <- df[df$strand == "+",]
df_minus <- df[df$strand == "-",]
#df <- df[with(df, order(chromosome, start, strand)),]


#####Load gene object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- as.vector(unique(df$chromosome))
tx_genes <- genes(txdb)


#####Define regions for TSSi as start/1000 (segments can't be wider than 1000 bp, important for computation time)
df_plus$region <- as.integer(df_plus$start / 1000)
df_minus$region <- as.integer(df_minus$start / 1000)
#make dataframe for TSSi
df_tssi_plus <- data.frame("chromosome" = df_plus$chromosome, "region" = df_plus$region, "start" = df_plus$start, "strand" = df_plus$strand, "counts" = df_plus$count)
df_tssi_minus <- data.frame("chromosome" = df_minus$chromosome, "region" = df_minus$region, "start" = df_minus$start, "strand" = df_minus$strand, "counts" = df_minus$count)
#drop NA rows, as segmentizeCounts reads them as region
df_tssi_plus <- df_tssi_plus[complete.cases(df_tssi_plus),]
df_tssi_minus <- df_tssi_minus[complete.cases(df_tssi_minus),]


#####TSSi workflow
###create segments
attach(df_tssi_plus)
xplus <- segmentizeCounts(counts=counts, start=start, chr=chromosome, region=region, strand=strand)
detach(df_tssi_plus)
attach(df_tssi_minus)
xminus <- segmentizeCounts(counts=counts, start=start, chr=chromosome, region=region, strand=strand)
detach(df_tssi_minus)

###normalize data in segments and split computation on several cores
yFit <- normalizeCounts(x, fit=TRUE, multicore = TRUE, mc.cores=12)
###identifyStartSites for each segments and split computation on several cores
z <- identifyStartSites(yFit, multicore = TRUE, mc.cores = 4)


#####export computed transcriptional start sites after annotation of regions
dfs <- tss(z)
tss_df <- bind_rows(dfs, .id ="identifier") 
###manipulate dataframe to contain information needed for transformation to GRanges (chromosome, start, end, strand )
tss_df$chromosome <- sapply(strsplit(tss_df$identifier, "_"), "[", 1)
tss_df$strand <- sapply(strsplit(tss_df$identifier, "_"), "[", 2)
tss_df$end <- tss_df$pos
tss_df <- rename(tss_df, c("pos" = "start"))
###create GRanges object and compute overlap to annotate genes
tss_ranges <- makeGRangesFromDataFrame(tss_df, keep.extra.columns = TRUE)
tss_df$overlaps=findOverlaps(tss_ranges, tx_genes, select = 'first')
tss_df$regions= with(tss_df, names(tx_genes)[overlaps])

###export data as csv
tss_export <- data.frame("geneID"=tss_df$regions, "Position"=tss_df$start, "BarcodeCount"=tss_df$reads, "ReadCount"=tss_df$reads, "PosFromAnno"="notdef", "Class"="notdef")

write.table(tss_export, file = "/home/stroemic/hiwi_16/analysis/TSSi/Fantom5/allFantomConsensusTss.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
