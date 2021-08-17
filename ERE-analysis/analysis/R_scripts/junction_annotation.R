library(tidyverse)
library(GenomicRanges)

import_directory='~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/External_packages/LIONS/mTEC-analysis'

junction_files = list.files(path=import_directory, pattern="*junctions.bed", 
                            full.names=TRUE, 
                            recursive=FALSE)
ensembl = useEnsembl(biomart = "genes", dataset = 'hsapiens_gene_ensembl')

transcript_annotation = getBM(mart = ensembl, attributes = c('ensembl_transcript_id',
                                                             'transcript_start',
                                                             'transcript_end',
                                                             'strand',
                                                             'chromosome_name',
                                                             'gene_biotype',
                                                             'external_gene_name',
                                                             'ensembl_gene_id'))

transcript_annotation = subset(transcript_annotation, chromosome_name %in% c(1:19, 'X', 'Y', 'MT')) %>%
  mutate(chromosome_name = glue::glue('chr{chromosome_name}')) %>%
  mutate(strand = case_when(strand == 1 ~ '+',
                            strand == -1 ~ '-'))

transcript_annotation = makeGRangesFromDataFrame(df = transcript_annotation, keep.extra.columns = T,
                                                 seqnames.field = 'chromosome_name', start.field = 'transcript_start',
                                                 end.field = 'transcript_end')

## Generate LTR annotation

TE_annotation = read.table(file = '~/Downloads/hg38_rmsk',
                           header = F, sep="\t",stringsAsFactors=FALSE, quote="")
TE_annotation = dplyr::rename(TE_annotation, 
                              chromosome = V6, 
                              start = V7, 
                              end = V8,
                              strand = V10,
                              name = V11,
                              class = V12,
                              family = V13) %>%
  dplyr::select(c('chromosome', 'start', 'end', 'strand', 'name', 'class', 'family')) %>%
  mutate(chromosome = case_when(chromosome == 'chrM' ~ 'chrMT',
                                T ~ chromosome))
TE_annotation = makeGRangesFromDataFrame(df = TE_annotation, 
                                         keep.extra.columns = T)

LTR_annotation = subset(TE_annotation, class != 'DNA')

counter = 0
for (file in junction_files){
  
  counter = counter + 1
  
  sample_name = unlist(stringr::str_split(file, pattern = '/'))[11]

  ## Import a list of splice junctions from TopHat output
  
  junctions = read.csv(file = junction_files[counter],
                       header = F,
                       sep = '\t') %>%
    rename(chromosome = V1,
           start = V2,
           end = V3,
           junction_ID = V4,
           score = V5,
           strand = V6) %>%
    dplyr::select(c('chromosome', 'start', 'end', 'junction_ID', 'score', 'strand')) %>%
    mutate(chromosome = case_when(chromosome == 'chrM' ~ 'chrMT',
                                  T ~ chromosome))
  
  junctions = makeGRangesFromDataFrame(junctions, keep.extra.columns = T)
  
  junctions_start = junctions
  end(junctions_start) = start(junctions)
  
  junctions_end = junctions
  start(junctions_end) = end(junctions)
  
  ## Find instances in which a splice junction starts in an LTR
  
  query = subset(junctions_start, strand == '+')
  subject = LTR_annotation
  
  overlaps_LTR_plus = findOverlaps(query = query, subject = subject)
  overlaps_LTR_plus = as.data.frame(overlaps_LTR_plus) %>%
    mutate(junction_ID = query[queryHits, ]$junction_ID) %>%
    mutate(LTR = subject[subjectHits, ]$name) %>%
    dplyr::select(c('junction_ID', 'LTR'))
  
  query = subset(junctions_end, strand == '-')
  subject = LTR_annotation
  
  overlaps_LTR_negative = findOverlaps(query = query, subject = subject)
  overlaps_LTR_negative = as.data.frame(overlaps_LTR_negative) %>%
    mutate(junction_ID = query[queryHits, ]$junction_ID) %>%
    mutate(LTR = subject[subjectHits, ]$name) %>%
    dplyr::select(c('junction_ID', 'LTR'))
  
  overlaps_LTR = rbind(overlaps_LTR_plus, overlaps_LTR_negative)
  
  ## Find instances in which a splice junction ends in a gene
  
  query = subset(junctions_end, strand == '+')
  subject = transcript_annotation
  
  overlaps_transcript_plus = findOverlaps(query = query, subject = subject)
  overlaps_transcript_plus = as.data.frame(overlaps_transcript_plus) %>%
    mutate(junction_ID = query[queryHits, ]$junction_ID) %>%
    mutate(ensembl_transcript_id = subject[subjectHits, ]$ensembl_transcript_id) %>%
    mutate(ensembl_gene_id = subject[subjectHits, ]$ensembl_gene_id) %>%
    dplyr::select(c('junction_ID', 'ensembl_transcript_id', 'ensembl_gene_id'))
  
  query = subset(junctions_start, strand == '-')
  subject = transcript_annotation
  
  overlaps_transcript_negative = findOverlaps(query = query, subject = subject)
  overlaps_transcript_negative = as.data.frame(overlaps_transcript_negative) %>%
    mutate(junction_ID = query[queryHits, ]$junction_ID) %>%
    mutate(ensembl_transcript_id = subject[subjectHits, ]$ensembl_transcript_id) %>%
    mutate(ensembl_gene_id = subject[subjectHits, ]$ensembl_gene_id) %>%
    dplyr::select(c('junction_ID', 'ensembl_transcript_id', 'ensembl_gene_id'))
  
  overlaps_transcript = rbind(overlaps_transcript_plus, overlaps_transcript_negative)
  
  ## Filter to only include junctions that also start in a TE
  
  overall_overlaps = subset(overlaps_LTR, junction_ID %in% overlaps_transcript$junction)
  
  transcript_name = vector()
  gene_name = vector()
  for (entry in 1:nrow(overall_overlaps)){
    
    transcripts = unique(subset(overlaps_transcript, junction_ID == overall_overlaps[entry, 'junction_ID'])$ensembl_transcript_id)
    transcript_name[entry] = stringr::str_c(transcripts,collapse = ';')
    
    genes = unique(subset(overlaps_transcript, junction_ID == overall_overlaps[entry, 'junction_ID'])$ensembl_gene_id)
    gene_name[entry] = stringr::str_c(genes,collapse = ';')
    
  }
  
  overall_overlaps$transcripts = transcript_name
  overall_overlaps$genes = gene_name
  
  ##
  
  junction_output = merge(x = overall_overlaps, y = as.data.frame(junctions), by = 'junction_ID') %>%
    mutate(sample = sample_name)
  
  if (counter == 1){
    
    aggregate_junction_output = junction_output
    
  }
  
  else{
    
    aggregate_junction_output = bind_rows(aggregate_junction_output, junction_output)
    
  }
  
}

int_only = aggregate_junction_output[grep(pattern = 'int', x = aggregate_junction_output$LTR), ]
int_only = int_only[order(int_only$genes), ]


################################################################################
## STAR junctions
################################################################################

ensembl = useEnsembl(biomart = "genes", dataset = 'mmusculus_gene_ensembl')

## Generate transcript annotation

transcript_annotation = getBM(mart = ensembl, attributes = c('ensembl_transcript_id',
                                                             'transcript_start',
                                                             'transcript_end',
                                                             'strand',
                                                             'chromosome_name',
                                                             'gene_biotype',
                                                             'external_gene_name',
                                                             'ensembl_gene_id'))

transcript_annotation = subset(transcript_annotation, chromosome_name %in% c(1:19, 'X', 'Y')) %>%
  mutate(chromosome_name = glue::glue('chr{chromosome_name}')) %>%
  mutate(strand = case_when(strand == 1 ~ '+',
                            strand == -1 ~ '-'))

transcript_annotation = makeGRangesFromDataFrame(df = transcript_annotation, keep.extra.columns = T,
                                                 seqnames.field = 'chromosome_name', start.field = 'transcript_start',
                                                 end.field = 'transcript_end')

## Generate TE annotation

TE_annotation = read.table(file = '~/tRF_targets/analysis/import/mm39_rmsk',
                           header = F, sep="\t",stringsAsFactors=FALSE, quote="")
TE_annotation = dplyr::rename(TE_annotation, 
                              chromosome = V6, 
                              start = V7, 
                              end = V8,
                              strand = V10,
                              name = V11,
                              class = V12,
                              family = V13) %>%
  dplyr::select(c('chromosome', 'start', 'end', 'strand', 'name', 'class', 'family'))
TE_annotation = makeGRangesFromDataFrame(df = TE_annotation, 
                                         keep.extra.columns = T)

LTR_annotation = subset(TE_annotation, class == 'LTR')

## Import the STAR juncions annotation

import_directory='/Users/mpeacey/tRF_targets/analysis/import/Macfarlan/STAR_junctions'

junction_files = list.files(path=import_directory, pattern="*SJ.out.tab", 
                            full.names=TRUE, 
                            recursive=FALSE)


canonical_chromosomes = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                          'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                          'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                          'chr19', 'chrX', 'chrY')
counter = 0
for (sample in junction_files){
  
  counter = counter + 1
  print(sample)
  print(counter)
  
  junctions = read.csv(file = sample,
                       header = F,
                       sep = '\t') %>%
    rename(chromosome = V1,
           start = V2,
           end = V3,
           strand = V4,
           intron_motif = V5,
           annotated= V6,
           unique_reads = V7,
           multi_reads = V8,
           alignment_overhand = V9) %>%
    mutate(strand = case_when(strand == 1 ~ '+',
                              strand == 2 ~ '-',
                              T ~ 'unknown')) %>%
    subset(chromosome %in% canonical_chromosomes) %>%
    subset(strand == '+' | strand == '-')
  
  junctions['junction_ID'] = rownames(junctions)
  
  junctions = makeGRangesFromDataFrame(junctions, keep.extra.columns = T)
  
  junctions_plus = subset(junctions, strand == '+')
  junctions_minus = subset(junctions, strand == '-')
  
  ## On the plus strand, identify instances in which a junction starts in an LTR 
  ## and ends in a gene
  
  junctions_plus_start = junctions_plus
  end(junctions_plus_start) = start(junctions_plus)
  
  LTR_overlaps = findOverlaps(query = junctions_plus_start, subject = LTR_annotation)
  LTR_overlaps = as.data.frame(LTR_overlaps) %>%
    mutate(junction_ID = junctions_plus_start[queryHits, ]$junction_ID) %>%
    mutate(LTR_name = LTR_annotation[subjectHits, ]$name) %>%
    dplyr::select(c('junction_ID', 'LTR_name'))
  
  junctions_plus_end = junctions_plus
  start(junctions_plus_end) = end(junctions_plus)
  
  transcript_overlaps = findOverlaps(query = junctions_plus_end, subject = transcript_annotation)
  transcript_overlaps = as.data.frame(transcript_overlaps) %>%
    mutate(junction_ID = junctions_plus_end[queryHits, ]$junction_ID) %>%
    mutate(ensembl_transcript_id = transcript_annotation[subjectHits, ]$ensembl_transcript_id) %>%
    mutate(ensembl_gene_id = transcript_annotation[subjectHits, ]$ensembl_gene_id) %>%
    dplyr::select(c('junction_ID', 'ensembl_transcript_id', 'ensembl_gene_id'))
  
  overall_overlaps = subset(LTR_overlaps, junction_ID %in% transcript_overlaps$junction_ID)
  
  transcript_name = vector()
  gene_name = vector()
  for (entry in 1:nrow(overall_overlaps)){
    
    transcripts = unique(subset(transcript_overlaps, junction_ID == overall_overlaps[entry, 'junction_ID'])$ensembl_transcript_id)
    transcript_name[entry] = stringr::str_c(transcripts,collapse = ';')
    
    genes = unique(subset(transcript_overlaps, junction_ID == overall_overlaps[entry, 'junction_ID'])$ensembl_gene_id)
    gene_name[entry] = stringr::str_c(genes,collapse = ';')
    
  }
  
  overall_overlaps$transcripts = transcript_name
  overall_overlaps$genes = gene_name
  
  df_plus = merge(x = overall_overlaps, y = junctions_plus, by = 'junction_ID') 
  
  ## On the minus strand, identify instances in which a junction starts in a gene 
  ## and ends an LTR
  
  junctions_minus_start = junctions_minus
  end(junctions_minus_start) = start(junctions_minus)
  
  transcript_overlaps = findOverlaps(query = junctions_minus_start, subject = transcript_annotation)
  transcript_overlaps = as.data.frame(transcript_overlaps) %>%
    mutate(junction_ID = junctions_minus_start[queryHits, ]$junction_ID) %>%
    mutate(ensembl_transcript_id = transcript_annotation[subjectHits, ]$ensembl_transcript_id) %>%
    mutate(ensembl_gene_id = transcript_annotation[subjectHits, ]$ensembl_gene_id) %>%
    dplyr::select(c('junction_ID', 'ensembl_transcript_id', 'ensembl_gene_id'))
  
  junctions_minus_end = junctions_minus
  start(junctions_minus_end) = end(junctions_minus)
  
  LTR_overlaps = findOverlaps(query = junctions_minus_end, subject = LTR_annotation)
  LTR_overlaps = as.data.frame(LTR_overlaps) %>%
    mutate(junction_ID = junctions_minus_end[queryHits, ]$junction_ID) %>%
    mutate(LTR_name = LTR_annotation[subjectHits, ]$name) %>%
    dplyr::select(c('junction_ID', 'LTR_name'))
  
  overall_overlaps = subset(LTR_overlaps, junction_ID %in% transcript_overlaps$junction_ID)
  
  transcript_name = vector()
  gene_name = vector()
  for (entry in 1:nrow(overall_overlaps)){
    
    transcripts = unique(subset(transcript_overlaps, junction_ID == overall_overlaps[entry, 'junction_ID'])$ensembl_transcript_id)
    transcript_name[entry] = stringr::str_c(transcripts,collapse = ';')
    
    genes = unique(subset(transcript_overlaps, junction_ID == overall_overlaps[entry, 'junction_ID'])$ensembl_gene_id)
    gene_name[entry] = stringr::str_c(genes,collapse = ';')
    
  }
  
  overall_overlaps$transcripts = transcript_name
  overall_overlaps$genes = gene_name
  
  df_minus = merge(x = overall_overlaps, y = junctions_minus, by = 'junction_ID') 
  
  df = bind_rows(df_plus, df_minus) %>%
    mutate(total_reads = unique_reads + multi_reads) %>%
    mutate(sample_name = sample)
  
  if (counter == 1){
    
    output = df
    
  }
  
  else{
    
    output = bind_rows(output, df)
    
  }
  
}

test = output[grep(pattern = 'int', x = output$LTR_name), ]

test = test[order(test$total_reads, decreasing = T), ]

