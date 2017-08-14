library(Gviz)

mTEC_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered.csv", header=TRUE, sep=",")
mTEC_pos_fil[mTEC_pos_fil$strand == "-",]$BarcodeCount <- (mTEC_pos_fil[mTEC_pos_fil$strand == "-",]$BarcodeCount)*-1
mTEC_pos_fil$strand <- "+"
mTEC_ranges <- makeGRangesFromDataFrame(mTEC_pos_fil, keep.extra.columns = TRUE)

plus <- mTEC_ranges[strand(mTEC_ranges) == '+']
minus <- mTEC_ranges[strand(mTEC_ranges) == '-']

mtrack <- DataTrack(mTEC_ranges, type="h")
mtrack_plus <- DataTrack(plus, type = "h")
mtrack_minus <- DataTrack(minus, type = "h")

chromosome = "chr11"

afrom = 18430000
ato = 18473859

ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chromosome)
axisTrack <- GenomeAxisTrack()


## Create genetrack from UCSC
# fill = "#8282d2" for colours 

genetrack=UcscTrack(track="RefSeq Genes", table="refGene", trackType="GeneRegionTrack",chromosome=chromosome, genome="hg19",
                    rstart="exonStarts", rends="exonEnds", gene="name", symbol="name2", transcript="name", strand="strand", 
                    stacking = "dense",
                    name="RefSeq Genes", feature="name2", showId=T, from = afrom, to = ato)



#plot_file_name = ""
# pdf(plot_file_name, width=12)
plotTracks(list(ideoTrack, axisTrack, genetrack, mtrack), from = afrom, to = ato,
           showBandId = TRUE, transcriptAnnotation = "symbol", #transcriptAnnotation= "symbol" for gene view
           groupAnnotation = "id")
# trash=dev.off()


biomTrack <- BiomartGeneRegionTrack(genome = "hg19", name = "ENSEMBL", symbol = "USHBP1", 
                                    transcriptAnnotation = "symbol", filters=list(hgnc_symbol="USHBP1"))

plotTracks(biomTrack)
