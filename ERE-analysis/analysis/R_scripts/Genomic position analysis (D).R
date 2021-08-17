library(DESeq2)
library(tidyverse)
library(GenomicRanges)
library(regioneR)
library(pheatmap)
library(RColorBrewer)
library(glue)

functions_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_functions/"
functions = c('extract_subset', 
              'differential_expression', 
              'process_DESeq2_results', 
              'make_GRanges', 
              'run_perm_test', 
              'calculate_odds_ratio_of_overlap',
              'save_heatmap_png')

for (i in functions){
  
  load(glue('{functions_directory}{i}'))
  
}

## Import required variables
results_df_local_ERE = readRDS(file = glue('{working_directory}/R_variables/results_df_local_ERE'))

#################################################################
# Generating GRanges objects
#################################################################

## TEs

GRanges_TE = make_GRanges(mode = 'TE', results_df = results_df_local_TE)

GRanges_TE_start = GRanges_TE
end(GRanges_TE_start) = start(GRanges_TE_start)

## Genes

gene_annotation = read.table(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/gencode.v38_gene_annotation_table.txt', header = 1)
gene_annotation = dplyr::select(gene_annotation, c('Geneid', 'Chromosome', 'Start', 'End', 'Strand', 'Class'))
gene_annotation = dplyr::rename(gene_annotation, chr = Chromosome, start = Start, end = End, strand = Strand)
gene_annotation$Geneid = gsub('\\..+$', '', gene_annotation$Geneid)
GRanges_gene = makeGRangesFromDataFrame(gene_annotation, keep.extra.columns = T)

saveRDS(GRanges_gene, file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_gene')

#GRanges_gene = make_GRanges(mode = 'gene', results_df = results_df_local_gene)

GRanges_gene_extended = GRanges_gene
end(GRanges_gene_extended) = end(GRanges_gene) + 1000
start(GRanges_gene_extended) = start(GRanges_gene) - 1000

saveRDS(GRanges_gene_extended, file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_gene_extended')

#GRanges_ERE = make_GRanges(mode = 'TE', results_df = results_df_local_ERE)

annotation = read.table(file = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt", header = 1)
annotation = separate(annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':')
annotation = separate(annotation, start.stop, into = c('start', 'end'), sep = '-')
annotation = dplyr::rename(annotation, locus = TE)

count_table_directory = glue::glue("{working_directory}/count_tables/TE_local")
list_of_EREs = extract_subset(read.table(glue::glue('{count_table_directory}/TE_local_hi_vs_lo.cntTable'),header=T,row.names=1), mode = 'ERE')
list_of_EREs$ID = row.names(list_of_EREs)
list_of_EREs = separate(list_of_EREs, col = 'ID', into = 'locus', sep = ':')$locus

ERE_annotation = subset(annotation, locus %in% list_of_EREs)

GRanges_ERE = makeGRangesFromDataFrame(ERE_annotation, keep.extra.columns = T)
saveRDS(GRanges_ERE, file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_ERE')

GRanges_ERE_start = GRanges_ERE
end(GRanges_ERE_start) = start(GRanges_ERE_start)

saveRDS(GRanges_ERE_start, file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_ERE_start')

GRanges_detected_ERE_start = GRanges_detected_EREs
end(GRanges_detected_ERE_start) = start(GRanges_detected_EREs)

# Genes of interest

AIRE_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/human_aire_dep_genes_san.csv')
FEZF2_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/human_fezf2_dep_genes.csv')
TRA_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/tra_genes.csv')
housekeeping_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/housekeeping_genes.csv')

up_genes = subset(results_df_local_gene, significant == T & log2FoldChange > 0)$Geneid
down_genes = subset(results_df_local_gene, significant == T & log2FoldChange < 0)$Geneid

# TEs of interest

up_EREs = subset(results_df_local_ERE, significant == T & log2FoldChange > 0)$locus
down_EREs = subset(results_df_local_ERE, significant == T & log2FoldChange < 0)$locus

#################################################################
# Gene differential expression (from Jason)
#################################################################

Sleuth_results_sig = read.csv('~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/Sleuth_results_sig.txt')

up_genes = subset(Sleuth_results_sig, b < 0)$ens_gene
down_genes = subset(Sleuth_results_sig, b > 0)$ens_gene

saveRDS(up_genes, file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/up_genes')
saveRDS(down_genes, file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/down_genes')

#################################################################
# ereMAPs
#################################################################

ereMAPs = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/ereMAPs_Larouche.csv', header = T)
ereMAPs$Start <- as.numeric(gsub(",","",ereMAPs$Start))
ereMAPs$End <- as.numeric(gsub(",","",ereMAPs$End))

GRanges_ereMAPs = makeGRangesFromDataFrame(ereMAPs, keep.extra.columns = T)

GRanges_detected_EREs = make_GRanges(results_df = results_df_local_ERE, mode = 'TE')

TE_ereMAP_overlaps = findOverlaps(query = GRanges_detected_EREs,
                                 subject = GRanges_ereMAPs)

TE_ereMAP_overlaps = data.frame(TE_ereMAP_overlaps)

validated = vector()
ereMAP_loci = vector()
for (entry in 1:nrow(TE_ereMAP_overlaps)){
  
  TE_hit = TE_ereMAP_overlaps[entry, ]$queryHits
  ereMAP_hit = TE_ereMAP_overlaps[entry, ]$subjectHits
  
  if (GRanges_detected_EREs[TE_hit, ]$gene == GRanges_ereMAPs[ereMAP_hit, ]$ERE.family){
    
    validated[entry] = T
    
  }
  
  else{
    
    validated[entry] = F

  }
  
  ereMAP_loci[entry] = GRanges_detected_EREs[TE_hit, ]$ID
  
}

results_df_local_ERE = mutate(results_df_local_ERE, ereMAP = case_when(ID %in% ereMAP_loci ~ T,
                                                                      !(ID %in% ereMAP_loci) ~ F))

output = generate_contingency(input = results_df_local_ERE,
                              condition_A = list('up_TEs' = 'significant == T & log2FoldChange > 0',
                                                 'down_TEs' = 'significant == T & log2FoldChange < 0',
                                                 'unchanged_TEs' = 'significant == F'),
                              condition_B = list('ereMAP' = 'ereMAP == T', 'not_ereMAP' = 'ereMAP == F'))

vcd::mosaic(~condition_A+condition_B, data = output, direction = c('v', 'h'), shade = T)


