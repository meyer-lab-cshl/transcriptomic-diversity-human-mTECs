library(GenomicRanges)
library(Gviz)
library(data.table)

gen <- "hg38"
chr <- "chr7"

afrom <- 150321976
ato <- 150323759

gene_track = AnnotationTrack(subset(GRanges_gene, Geneid == 'ENSG00000127399' | Geneid == 'ENSG00000188707'), name = 'genes')
TE_track = AnnotationTrack(subset(GRanges_ERE, locus == 'MER41B_dup1196'), name = 'TEs')
gtrack <- GenomeAxisTrack()
itrack = IdeogramTrack(genome = gen, chromosome = chr)

## Read track

files = list.files(path="~/Desktop/IGB_test", 
                   pattern="*new.bedgraph", 
                   full.names=TRUE, recursive=FALSE)
counter = 1
tracks = vector(mode = 'list')
for (file in files){
  
  bed <- fread(file, col.names = c('chromosome', 'start', 'end', 'value'))
  bed <- bed[chromosome == chr]
  
  if (grepl('hi', file) == T){color = '#4c72b0ff'}
  else{color = '#dd8452ff'}
  
  tracks[counter] = DataTrack(range = bed, type = "h", genome = gen, name = "Seq. Depth",
                    col = color, ylim = c(0,150))
  counter = counter + 1
  
}

## Plot

track_list = list(gtrack, tracks[[1]], tracks[[3]], tracks[[5]], tracks[[2]],
                  tracks[[4]], tracks[[6]], gene_track, TE_track)

plotTracks(track_list, from = afrom, to = ato)

read.table("~/Desktop/IGB_test/pt214_mTEC-hi_new.bedgraph.gz",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

exonskip_gtex <- import.bed('~/Desktop/IGB_test/pt214_mTEC-hi_new.bedgraph.gz')

plotTracks(list(gtrack, TE_track, gene_track))