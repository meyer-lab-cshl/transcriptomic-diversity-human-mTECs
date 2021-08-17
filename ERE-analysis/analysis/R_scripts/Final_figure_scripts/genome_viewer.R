library(GenomicRanges)
library(Gviz)
library(data.table)

GRanges_gene = readRDS(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_gene')
GRanges_ERE = readRDS(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_ERE')

gen <- "hg38"
chr <- "chr7"

TE_track = AnnotationTrack(subset(GRanges_ERE, locus == 'MER41B_dup1196'), name = 'TEs')
gtrack <- GenomeAxisTrack()
itrack = IdeogramTrack(genome = gen, chromosome = chr)

from <- 150322000
to <- 150338226

gtrack <- UcscTrack(genome = gen, chromosome = chr, 
                                  track = "All GENCODE V38", from = from, to = to,
                                  trackType = "GeneRegionTrack", 
                                  rstarts = "exonStarts", rends = "exonEnds", 
                                  gene = "name", symbol = "name", 
                                  transcript = "name", strand = "strand", 
                                  fill = "#8282d2", name = "UCSC Genes")

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
  
  tracks[counter] = DataTrack(range = bed, type = "h", genome = gen, name = file,
                    col = color, ylim = c(0,200))
  counter = counter + 1
  
}

## Plot

track_list = list(tracks[[1]], tracks[[3]], tracks[[5]], tracks[[2]],
                  tracks[[4]], tracks[[6]], gtrack, TE_track)

plotTracks(track_list, from = afrom, to = ato)

read.table("~/Desktop/IGB_test/pt214_mTEC-hi_new.bedgraph.gz",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

exonskip_gtex <- import.bed('~/Desktop/IGB_test/pt214_mTEC-hi_new.bedgraph.gz')

plotTracks(list(gtrack, TE_track, gene_track))