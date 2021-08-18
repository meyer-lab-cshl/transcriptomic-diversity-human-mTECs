library(GenomicRanges)
library(Gviz)
library(data.table)
library(glue)

working_directory = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis'

##Import required variables
GRanges_ERE = readRDS(file = glue('{working_directory}/R_variables/GRanges_ERE'))

##Track settings
gen = "hg38"
chr = "chr7"
from = 150322000
to = 150324000

## Ideogram track
itrack = IdeogramTrack(genome = gen, chromosome = chr, fontsize=7)

## Genome axis track
gtrack = GenomeAxisTrack()

## Transcript track
gene_track = UcscTrack(genome = gen, chromosome = chr, 
                        track = "All GENCODE V38", from = from, to = to,
                        trackType = "GeneRegionTrack", 
                        rstarts = "exonStarts", rends = "exonEnds", 
                        gene = "name", symbol = "name", 
                        transcript = "name", strand = "strand", 
                        fill = "#CCEBC5", name = "GENCODE transcripts",
                        showId = TRUE, col = 'black')

counter = 0
for (entry in symbol(gene_track)){
  
  counter = counter + 1
  symbol(gene_track)[counter] = glue::glue('LRRC61 ({entry})')
  
}

## TE track
TE_track = AnnotationTrack(subset(GRanges_ERE, locus == 'MER41B_dup1196'), 
                           name = 'EREs',
                           id = c('MER41B'), 
                           groupAnnotation='id',
                           showId=T,
                           fill = "#FBB4AE",
                           col = 'black',
                           col.id="black",
                           fontcolor.item='black',
                           showTitle = F)

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
  
  tracks[counter] = DataTrack(range = bed, type = "histogram", genome = gen, name = file,
                    fill = color, ylim = c(0,1), col.histogram = 'black', alpha = 0.9)
  counter = counter + 1
  
}

## Plot
track_list = list(itrack, gtrack, gene_track, TE_track, 
                  tracks[[1]], tracks[[3]], tracks[[5]], tracks[[2]],
                  tracks[[4]], tracks[[6]])

plotTracks(track_list, from = from, to = to, showTitle = FALSE, showAxis=F)