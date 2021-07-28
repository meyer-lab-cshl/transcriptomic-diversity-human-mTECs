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

#################################################################
# Generating GRanges objects
#################################################################

## TEs

GRanges_TE = make_GRanges(mode = 'TE', results_df = results_df_local_TE)

GRanges_TE_start = GRanges_TE
end(GRanges_TE_start) = start(GRanges_TE_start)

## Genes

annotation = read.table(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/gencode.v38_gene_annotation_table.txt', header = 1)
annotation = select(annotation, c('Geneid', 'Chromosome', 'Start', 'End', 'Strand', 'Class'))
annotation = dplyr::rename(annotation, chr = Chromosome, start = Start, end = End, strand = Strand)
annotation$Geneid = gsub('\\..+$', '', annotation$Geneid)
GRanges_gene = makeGRangesFromDataFrame(annotation, keep.extra.columns = T)

#GRanges_gene = make_GRanges(mode = 'gene', results_df = results_df_local_gene)

GRanges_gene_extended = GRanges_gene
end(GRanges_gene_extended) = end(GRanges_gene) + 1000
start(GRanges_gene_extended) = start(GRanges_gene) - 1000

saveRDS(GRanges_gene_extended, file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_gene_extended')

#GRanges_gene_upstream = GRanges_gene
#end(GRanges_gene_upstream) = start(GRanges_gene)
#start(GRanges_gene_upstream) = start(GRanges_gene) - 5000

## EREs

#GRanges_ERE = make_GRanges(mode = 'TE', results_df = results_df_local_ERE)

annotation = read.table(file = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt", header = 1)
annotation = separate(annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':')
annotation = separate(annotation, start.stop, into = c('start', 'end'), sep = '-')
annotation = dplyr::rename(annotation, locus = TE)

list_of_EREs = extract_subset(read.table(glue::glue('{count_table_directory}TE_local_hi_vs_lo.cntTable'),header=T,row.names=1), mode = 'ERE')
list_of_EREs$ID = row.names(list_of_EREs)
list_of_EREs = separate(list_of_EREs, col = 'ID', into = 'locus', sep = ':')$locus

ERE_annotation = subset(annotation, locus %in% list_of_EREs)

GRanges_ERE = makeGRangesFromDataFrame(ERE_annotation, keep.extra.columns = T)

GRanges_ERE_start = GRanges_ERE
end(GRanges_ERE_start) = start(GRanges_ERE_start)

saveRDS(GRanges_ERE_start, file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/GRanges_ERE_start')

GRanges_detected_ERE_start = GRanges_detected_EREs
end(GRanges_detected_ERE_start) = start(GRanges_detected_EREs)

# Genes of interest

AIRE_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/human_aire_dep_genes_san.csv')$ensembl_gene_id
FEZF2_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/human_fezf2_dep_genes.csv')$ensembl_gene_id
TRA_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/tra_genes.csv')$ensembl_gene_id

housekeeping_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/housekeeping_genes.csv')$ensembl_gene_id
housekeeping_genes = housekeeping_genes[!(housekeeping_genes %in% AIRE_genes)]
housekeeping_genes = housekeeping_genes[!(housekeeping_genes %in% FEZF2_genes)]

#ensembl = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
#housekeeping_genes = getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", values = housekeeping_genes, mart= ensembl)$ensembl_gene_id

up_genes = subset(results_df_local_gene, significant == T & log2FoldChange > 0)$Geneid
down_genes = subset(results_df_local_gene, significant == T & log2FoldChange < 0)$Geneid

# TEs of interest

up_EREs = subset(results_df_local_ERE, significant == T & log2FoldChange > 0)$locus
down_EREs = subset(results_df_local_ERE, significant == T & log2FoldChange < 0)$locus

###############
## TEeffectR ##
###############

library(stringr)
library(biomaRt)
library(biomartr)
library(dplyr)
library(Rsamtools)
library(edgeR)
library(rlist)
library(limma)
library(TEffectR)
library(tidyr)
library(GenomicRanges)

## Counts

gene.counts = extract_subset(mode = 'gene', input = data)

TE.counts = extract_subset(mode = 'TE', input = data) 
TE.counts = cbind(ID = rownames(TE.counts), TE.counts)
TE.counts = separate(data = TE.counts, col = 'ID', into = c('locus', 'gene', 'family', 'class'), sep = ':')
TE.counts = cbind(ID = rownames(TE.counts), TE.counts)
TE.counts = subset(TE.counts, class == 'LTR')

## Annotations

gene.annotation = read.table(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/gencode.v38_gene_annotation_table.txt', header = 1) %>%
  select(c('Chromosome', 'Start', 'End', 'Strand', 'Geneid', 'GeneSymbol')) %>%
  dplyr::rename(chr = Chromosome, start = Start, end = End, strand = Strand, geneID = Geneid, geneName = GeneSymbol)

repeat.annotation = read.table(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt', header = 1) %>%
  separate(chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':') %>%
  separate(start.stop, into = c('start', 'end'), sep = '-') %>%
  dplyr::rename(locus = TE)

#repeatmasker.annotation <- TEffectR::rm_format(filepath = "~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/hg38.fa.out")

## Make GRanges objects

TE.counts.GRanges = merge(TE.counts, repeat.annotation, by = 'locus') %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
TE.counts.GRanges.start = TE.counts.GRanges
end(TE.counts.GRanges.start) = start(TE.counts.GRanges)
start(TE.counts.GRanges.start) = start(TE.counts.GRanges)

gene.annotation.GRanges = makeGRangesFromDataFrame(gene.annotation, keep.extra.columns = T)
gene.annotation.GRanges.upstream = gene.annotation.GRanges
end(gene.annotation.GRanges.upstream) = start(gene.annotation.GRanges)
start(gene.annotation.GRanges.upstream) = start(gene.annotation.GRanges) - 5000

## Build sum.repeat.counts

query = gene.annotation.GRanges.upstream
subject = TE.counts.GRanges.start

hits = as.data.frame(findOverlaps(query, subject))

sum.repeat.counts = hits %>%
  mutate(queryHits = gene.annotation[queryHits, 'geneName']) %>%
  dplyr::rename(geneName = queryHits) %>%
  mutate(subjectHits = TE.counts[subjectHits, 2]) %>%
  dplyr::rename(locus = subjectHits) %>%
  merge(TE.counts, by = 'locus') %>%
  select(-c('ID', 'gene', 'family')) %>%
  dplyr::rename(repeatName = locus, repeatClass = class)

## RUn TEffectR

# read your gene annotation file
#gene.annotation<-read.table("~/Desktop/thymus-epitope-mapping/ERE-analysis/sampleInputs/gene.annotation.tsv", header= T, stringsAsFactors = F)

# read your gene expression file
#gene.counts<-read.table("~/Desktop/thymus-epitope-mapping/ERE-analysis/sampleInputs/gene.counts.tsv", header= T, row.names=1, stringsAsFactors = F)

# read your summarised repeat annotation file
#sum.repeat.counts<-read.table("~/Desktop/thymus-epitope-mapping/ERE-analysis/sampleInputs/sum.repeat.counts.tsv", header= T, stringsAsFactors = F)

covariates <- NULL
prefix = "LTRs"

lm = TEffectR::apply_lm(gene.annotation = gene.annotation,
                       gene.counts = gene.counts,
                       repeat.counts = sum.repeat.counts,
                       covariates = covariates,
                       prefix = prefix)

lm_results = read.table("~/Desktop/thymus-epitope-mapping/ERE-analysis/LTRs -lm-results.tsv", header= T, sep="\t") %>%
  mutate(significant = case_when(model.p.value < 0.05 ~ T,
                                 model.p.value >= 0.05 ~ F))


#  mutate(adjusted.p.value = p.adjust(model.p.value, method = 'BH', n = nrow(lm_results))) %>%
#  mutate(significant = case_when(adjusted.p.value < 0.05 ~ T,
#                                 adjusted.p.value >= 0.05 ~ F))

volcano_plot = ggplot() +
  geom_point(data = lm_results, aes(x = r.squared, y = -log10(model.p.value)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(lm_results, significant == TRUE), aes(x = r.squared, y = -log10(model.p.value), fill = significant), size = 2.5, alpha = 1, shape = 21, stroke = 0) +
  geom_point(data = subset(lm_results, significant == FALSE), aes(x = r.squared, y = -log10(model.p.value)), size = 1, alpha = 1, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('R squared')) +
  ylab(expression('-log'[10]*'(p-value)')) +
  ggrepel::geom_label_repel(data = subset(lm_results, significant == TRUE), aes(x = r.squared, y = -log10(model.p.value), label = GeneName))

volcano_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank(),
                                  legend.text = element_text(size = 15),
                                  legend.title = element_text(size = 18),
                                  legend.position = c(0.2, 0.93),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/TEeffectR_LTRs.png", 
       width = 7, height = 5.25, units = "in")

plot_regression = function(regression, counts, gene){
  
  counts = cbind(ID = counts$X, counts) %>%
    separate(col = 'ID', into = c('associated_gene', 'repeat_class', 'repeat_name'), sep = ':')

  x_values = subset(counts, X == gene) %>%
    select(-c('repeat_class', 'repeat_name', 'associated_gene')) %>%
    pivot_longer(cols = -1, values_to = 'gene_expression')
  
  TE = subset(regression, GeneName == gene)$RepeatName 
  
  y_values = subset(counts, repeat_name == TE) %>%
    select(-c('repeat_class', 'X', 'associated_gene')) %>%
    pivot_longer(cols = -1, values_to = 'TE_expression')
  
  output = merge(x_values, y_values, by = 'name') %>%
    separate(col = 'name', into = c('patient', 'mTEC'), sep = '_')

  return(output)
  
}

lm_cpm = read.table("~/Desktop/thymus-epitope-mapping/ERE-analysis/LTRs -cpm-values.tsv", header= T, sep="\t") 
lm_results = read.table("~/Desktop/thymus-epitope-mapping/ERE-analysis/LTRs -lm-results.tsv", header= T, sep="\t")

regression = lm_results
counts = lm_cpm
gene = 'DHDH'

correlation_df = plot_regression(regression = regression, counts = counts, gene = gene)

correlation_plot = ggplot(data = correlation_df, aes(x = gene_expression, y = TE_expression)) +
  geom_point(aes(fill = mTEC), shape = 21, size = 4) +
  stat_smooth(method = "lm", col = "red") +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff')) +
  xlab(expression('Gene Expression (log'[2]*'(CPM))')) +
  ylab(expression('TE Expression (log'[2]*'(CPM))'))

correlation_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                        plot.subtitle = element_text(size = 14),
                                        axis.text.x = element_text(size = 14),
                                        axis.text.y = element_text(size = 14),
                                        axis.title = element_text(size = 14),
                                        axis.line = element_line(size = 0.8),
                                        panel.border = element_blank(),
                                        legend.text = element_text(size = 15),
                                        legend.title = element_text(size = 18),
                                        legend.position = c(0.2, 0.93),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/TEeffectR_TBC1D15.png", 
       width = 7, height = 5.25, units = "in")



################
## Odds ratio ##
################

calculate_odds_ratio_of_overlap = function(gene_sets,
                                           TE_set,
                                           TE_coordinates = GRanges_TE, 
                                           gene_coordinates = GRanges_gene){
  
  ## 'gene sets' is a list in which each entry is a list of ENSEMBL gene IDs.
  ##
  ## 'TE_coordinates' is a GRanges object containing all detected TEs  and
  ## the output of differential expression analysis via TElocal.
  ##
  ## 'gene_coordinates' is a GRanges object containing all detected genes and
  ## and the output of differential expression analysis via TElocal.
  
  subset = names(gene_sets)
  p_value = vector()
  odds_ratio = vector()
  lower_interval = vector()
  upper_interval = vector()
  
  if (TE_set == 'up'){
    
    fold_change_condition = 'log2FoldChange > 0'
    
  }
  
  if (TE_set == 'down'){
    
    fold_change_condition = 'log2FoldChange < 0'
    
  }
  
  for (i in 1:length(gene_sets)){
    
    gene_set = gene_sets[[i]]
    
    column_1 = c(length(findOverlaps(subset(TE_coordinates, significant == T & eval(parse(text = fold_change_condition))), subset(gene_coordinates, Geneid %in% gene_set))), 
                 length(findOverlaps(subset(TE_coordinates, !(significant == T & eval(parse(text = fold_change_condition)))), subset(gene_coordinates, Geneid %in% gene_set))))
    
    column_2 = c(length(findOverlaps(subset(TE_coordinates, significant == T & eval(parse(text = fold_change_condition))), subset(gene_coordinates, !(Geneid %in% gene_set)))), 
                 length(findOverlaps(subset(TE_coordinates, !(significant == T & eval(parse(text = fold_change_condition)))), subset(gene_coordinates, !(Geneid %in% gene_set)))))
    
    contingency = data.frame('gene' = column_1, 'not gene' = column_2, row.names = c('TE', 'not TE'))
    
    p_value[i] = fisher.test(x = contingency)[[1]]
    odds_ratio[i] = fisher.test(x = contingency)[[3]]
    lower_interval[i] = fisher.test(x = contingency)$conf.int[1]
    upper_interval[i] = fisher.test(x = contingency)$conf.int[2]
    
  }
  
  p_value = p.adjust(p_value, method = 'bonferroni')
  
  output = data.frame(subset = subset, 
                      p_value = p_value, 
                      odds_ratio = odds_ratio, 
                      lower_interval = lower_interval, 
                      upper_interval = upper_interval) %>%
    mutate(significant = case_when(p_value < 0.001 ~ '***', 
                                   p_value < 0.01 ~ '**',
                                   p_value < 0.05 ~ '*',
                                   p_value >= 0.05 ~ '')) %>%
    mutate(subset = forcats::fct_reorder(subset, odds_ratio))
  
  return(output)
  
}

gene_sets = list('AIRE' = AIRE_genes, 
                 'FEZF2' = FEZF2_genes, 
                 'TRA' = TRA_genes,
                 'Housekeeping' = housekeeping_genes)

#gene_sets = list('up_genes' = subset(GRanges_gene, significant == T & log2FoldChange > 0)$Geneid,
#                 'down_genes' = subset(GRanges_gene, significant == T & log2FoldChange < 0)$Geneid)

output = calculate_odds_ratio_of_overlap(gene_sets = gene_sets,
                                         TE_set = 'down',
                                         TE_coordinates = GRanges_ERE_start,
                                         gene_coordinates = GRanges_gene_extended)

odds_ratio_plot = ggplot(data = output, aes(x = subset, y = odds_ratio)) + 
  geom_point() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_errorbar(aes(ymin = lower_interval, ymax = upper_interval), width = 0) +
  xlab('') +
  ylab('Odds ratio') +
  geom_text(aes(label = significant), nudge_y = 0.3, size = 6) +
  scale_y_continuous(trans='log10')

odds_ratio_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                     plot.subtitle = element_text(size = 14),
                                     axis.text.x = element_text(size = 14),
                                     axis.text.y = element_text(size = 14),
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_text(size = 15, margin = margin(r = 7.5)),
                                     axis.line = element_line(size = 0.8),
                                     panel.border = element_blank(),
                                     legend.text = element_text(size = 15),
                                     legend.title = element_text(size = 18),
                                     legend.position = c(0.2, 0.93),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/odds-ratio_overlap_genes-of-interest_down.png", 
       width = 5, height = 5, units = "in")


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


