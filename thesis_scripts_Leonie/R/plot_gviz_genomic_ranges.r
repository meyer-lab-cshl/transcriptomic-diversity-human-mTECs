library(Gviz)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

## functions
preprocess_df <- function(x) {
  x[x$strand == "-",]$BarcodeCount <- (x[x$strand == "-",]$BarcodeCount)*-1
  x$strand <- "+"
  return(x)
}

### read in 
mTEC_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_mTECs.filtered.csv", header=TRUE, sep=",")
tissues_pos_fil <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.tissues.filtered.csv", header=TRUE, sep=",")



##### BioMart object
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="grch37.ensembl.org", 
                   dataset="hsapiens_gene_ensembl")

##### txdb gene object
my_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)<- my_chr


#### pre
mTEC_pos_fil <- preprocess_df(mTEC_pos_fil)
mTEC_ranges <- makeGRangesFromDataFrame(mTEC_pos_fil, keep.extra.columns = TRUE)
tissues_pos_fil <- preprocess_df(tissues_pos_fil)
tissues_ranges <- makeGRangesFromDataFrame(tissues_pos_fil, keep.extra.columns = TRUE)


## list fo genes 
genes <- c("CTRB2", "CYP3A4","USHBP1","STAR","LDHC","CELA2B","RARRES1","CHIT1","RHCG","CD46","TG")
chromosomes <- list("CTRB2" = "chr16", "CYP3A4"="chr7","USHBP1"="chr19","STAR"="chr8","LDHC"="chr11",
                    "CELA2B"="chr1","RARRES1"="chr3","CHIT1"="chr1","RHCG"="chr15","CD46"="chr1","TG"="chr8")

entrezID <- list("CTRB2" ="440387", "CYP3A4"="1576","USHBP1"="83878","STAR"="6770","LDHC"="3948",
                 "CELA2B"="51032","RARRES1"="5918","CHIT1"="1118","RHCG"="51458","CD46"="4179","TG"="7038")
                    
starts <- list("CTRB2" = 75237550, "CYP3A4"=99351650,"USHBP1"=17359350,"STAR"=37998880,"LDHC"=18430380,
               "CELA2B"=15800650,"RARRES1"=158412450,"CHIT1"=203184000,"RHCG"=90012230,"CD46"=207921900,"TG"=133862210)
ends <- list("CTRB2" = 75241200, "CYP3A4"=99382900,"USHBP1"=17376000,"STAR"=38008780,"LDHC"=18473310,
             "CELA2B"=15819850,"RARRES1"=158451900,"CHIT1"=203199350,"RHCG"=90040450,"CD46"=207969400,"TG"=134150850)



#First Exon
exon = exonsBy(txdb, by='gene')
oe <- order(as.numeric(names(exon)))
exon <- exon[oe]
#Other Exons
other_exon = lapply(exon, function(x) if( unique(strand(x)) == '+') x[-1] else x[-length(x)])
other_exon = do.call(GRangesList, other_exon)
other_exon = unlist(other_exon)

#TSS10 & 100 other exon
TSS10_other_exon <- promoters(other_exon, upstream = 10, downstream = 10)
TSS100_other_exon <- promoters(other_exon, upstream = 100, downstream = 100)
strand(TSS100_other_exon) = "*"
mcols(TSS100_other_exon) = 1



for(i in genes[2]){
  chromosome <- chromosomes[[i]]
  afrom <- starts[[i]]
  ato <- ends[[i]]
  entrez <- entrezID[[i]]
  
  
  m_max <- max(mTEC_pos_fil[mTEC_pos_fil$chrom == chromosome & mTEC_pos_fil$start > afrom &  mTEC_pos_fil$start < ato ,]$BarcodeCount)
  if(m_max < 0) {m_max = 0}
  m_min <- min(mTEC_pos_fil[mTEC_pos_fil$chrom == chromosome & mTEC_pos_fil$start > afrom &  mTEC_pos_fil$start < ato ,]$BarcodeCount)
  if(m_min > 0) {m_min = 0}
  mtrack <- DataTrack(mTEC_ranges, type="h", baseline=0, name="mTECs", 
                      col.baseline="grey", ylim=c(m_min-(0.1*m_min),m_max+(0.1*m_max)))
  
  t_max <- max(tissues_pos_fil[tissues_pos_fil$chrom == chromosome & tissues_pos_fil$start > afrom &  tissues_pos_fil$start < ato ,]$BarcodeCount)
  if(t_max < 0) {t_max = 0}
  t_min <- min(tissues_pos_fil[tissues_pos_fil$chrom == chromosome & tissues_pos_fil$start > afrom &  tissues_pos_fil$start < ato ,]$BarcodeCount)
  if(t_min > 0) {t_min = 0}
  ttrack <- DataTrack(tissues_ranges, type="h", baseline=0, name="tissues wo thymus", 
                      col.baseline="grey", ylim=c(t_min-(0.1*t_min),t_max+(0.1*t_max)))
  
  ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chromosome)
  axisTrack <- GenomeAxisTrack()
  
  if(i==i){
    TSS100 <- TSS100_other_exon[names(TSS100_other_exon)==entrez]
    TSS100 <- reduce(TSS100)
    TSS100_track <- AnnotationTrack(TSS100, name="Exon TSS100", col="#800026", fill="#800026")
    
    ## Create genetrack from UCSC
    # fill = "#8282d2" for colours 
    for(x in c("dense","squish")){
      genetrack=UcscTrack(track="RefSeq Genes", table="refGene", trackType="GeneRegionTrack",chromosome=chromosome, genome="hg19",
                          stacking = x,
                          rstart="exonStarts", rends="exonEnds", gene="name", symbol="name2", transcript="name", strand="strand", 
                          name="RefSeq Genes", feature="name2", geneSymbol=T, from = afrom, to = ato)
      
      
      
      plot_file_name = paste("/home/stroemic/hiwi_16/analysis/shifted_TSS/plots/gviz/",i,"_gviz_range_",x,"_oe.pdf",sep="")
      pdf(plot_file_name, width=12)
      plotTracks(list(ideoTrack, axisTrack, genetrack, TSS100_track,ttrack, mtrack), from = afrom, to = ato,
                 transcriptAnnotation = "symbol")#transcriptAnnotation= "symbol" for gene view
      
      dev.off()
    }
  } else {
    ## Create genetrack from UCSC
    # fill = "#8282d2" for colours 
    for(x in c("dense","squish")){
      genetrack=UcscTrack(track="RefSeq Genes", table="refGene", trackType="GeneRegionTrack",chromosome=chromosome, genome="hg19",
                      stacking = x,
                      rstart="exonStarts", rends="exonEnds", gene="name", symbol="name2", transcript="name", strand="strand", 
                      name="RefSeq Genes", feature="name2", geneSymbol=T, from = afrom, to = ato)
  
  
  
      plot_file_name = paste("/home/stroemic/hiwi_16/analysis/shifted_TSS/plots/gviz/",i,"_gviz_range_",x,"_.pdf",sep="")
      pdf(plot_file_name, width=12)
      plotTracks(list(ideoTrack, axisTrack, genetrack, ttrack, mtrack), from = afrom, to = ato,
            transcriptAnnotation = "symbol")#transcriptAnnotation= "symbol" for gene view
      
      dev.off()
    }
  }
  
  
  
  
  
}




