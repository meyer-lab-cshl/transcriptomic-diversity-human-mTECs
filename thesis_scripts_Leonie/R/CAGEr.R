library(dplyr)
library(CAGEr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("/home/stroemic/hiwi_16/data/external/Fantom5/CTSS_BED/tissues/")


#####Load gene object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")
seqlevels(txdb)<- my_chr
tx_genes <- genes(txdb)

## read in file list, store names of tissues in list, append all_mTECs sample
fantoms <- list.files(pattern="*.bed$", full.names = TRUE)
samples <- lapply(fantoms, function(x) gsub("./(.*).CNhs.*", "\\1", x ))
fantoms <- append(fantoms, "/home/stroemic/hiwi_16/data/Summary_counts/bed/all_mTECs.bed")
samples <- append(samples, "all_mTECs")


###CAGEr workflow
###make CAGE dataset
myCAGEset <- new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg19",
                 inputFiles = fantoms, inputFilesType = "bed",
                 sampleLabels = unlist(samples))
##read in TSS positions and save because time intensive
getCTSS(myCAGEset)
saveRDS(object = myCAGEset, file = "/home/stroemic/hiwi_16/analysis/CAGEr/myCAGEset.rds")
## calculate TSS 
ctss_calc <- CTSStagCount(myCAGEset)
## normalisation --> here it failed for the big dataset, even in the simpleTpm mode
normalizeTagCount(myCAGEset, method = "simpleTpm")

## cluster TSS 
clusterCTSS(object = myCAGEset, threshold = 2, thresholdIsTpm = TRUE, nrPassThreshold = 1,  
            method = "distclu", maxDist = 12, removeSingletons = TRUE, keepSingletonsAbove = 10,
            useMulticore = TRUE, nrCores = 10)

cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters", useMulticore = TRUE, nrCores = 10)
quantilePositions(myCAGEset,clusters = "tagClusters", qLow=0.1, qUp=0.9, useMulticore = TRUE, nrCores = 10)


######################################################
# setwd("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/")
# ##read in filtered data files
# files <- list.files(pattern = "*ctss.filtered.csv$")
# files <- append(files, "all_mTECs.filtered.csv")
# 
# 
# c = 1
# dats <- list()
# for (i in files) {
#   dat <- read.table(i, sep = ",", header=TRUE)
#   if (grepl("mTEC", i)) {
#     tissue <- "all_mTECs"
#   } else {
#     tissue <- gsub("^(.*).CNhs.*", "\\1", i )
#   }
#   
#   print(tissue)
#   #dat <- dat[!grepl("intergenic", dat$geneID),]
#   dat_out <- data.frame("chr"=dat$chrom, "pos"=dat$start, "strand"=dat$strand,  "x"=as.integer(dat$BarcodeCount))
#   names(dat_out)[names(dat_out) == "x"] <- tissue
#   dats[[c]] <- dat_out
#   c = c+1
# }
# 
# 
# ctss <- Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by=c("chr","pos","strand")), dats)
# ctss[is.na(ctss)] <- as.integer(0)
# ctss$strand <- as.character(ctss$strand)
# 
# myCAGEset <- as(ctss, "CAGEset")
# 
# ctss_calc <- CTSStagCount(myCAGEset)
# 
# clusterCTSS(object = myCAGEset, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1, method="distclu", maxDist=20, removeSingletons=TRUE, keepSingletonsAbove=5)
