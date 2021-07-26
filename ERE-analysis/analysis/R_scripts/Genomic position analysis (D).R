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
# DESeq2
#################################################################

count_table_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/"
data = read.table(glue::glue('{count_table_directory}TE_local_hi_vs_lo.cntTable'),header=T,row.names=1)

dds_local = differential_expression(data, design=~patient+tissue)

results_local = results(dds_local, 
                        independentFiltering = F,
                        contrast = c('tissue', 'hi', 'lo'))

results_local_gene = extract_subset(mode = 'gene', input = results_local)
results_local_TE = extract_subset(mode = 'TE', input = results_local)
results_local_ERE = extract_subset(mode = 'ERE', input = results_local)

results_df_local_gene = process_DESeq2_results(results = results_local_gene, mode = 'Gene') 

keep_list = c('LTR', 'LINE', 'SINE', 'Satellite', 'DNA')

results_df_local_TE = process_DESeq2_results(results = results_local_TE, mode = 'TE_local') %>%
  mutate(ID = sub("\\?", "", ID)) %>%
  mutate(class = sub("\\?", "", class)) %>%
  mutate(class = case_when(class %in% keep_list ~ class,
                           !(class %in% keep_list) ~ 'Other'))

results_df_local_ERE = process_DESeq2_results(results = results_local_ERE, mode = 'TE_local') %>%
  mutate(ID = sub("\\?", "", ID)) %>%
  mutate(class = sub("\\?", "", class))

results_df_local_gene_up = filter(results_df_local_gene, (significant == T) & (log2FoldChange > 0))
results_df_local_gene_unchanged = filter(results_df_local_gene, significant == F)
results_df_local_gene_down = filter(results_df_local_gene, (significant == T) & (log2FoldChange < 0))
results_df_local_gene_sigdiff = filter(results_df_local_gene, significant == T)

results_df_local_TE_up = filter(results_df_local_TE, (significant == T) & (log2FoldChange > 0))
results_df_local_TE_unchanged = filter(results_df_local_TE, significant == F)
results_df_local_TE_down = filter(results_df_local_TE, (significant == T) & (log2FoldChange < 0))
results_df_local_TE_sigdiff = filter(results_df_local_TE, significant == T)

#################################################################
# Generating GRanges objects
#################################################################

## TEs

GRanges_TE = make_GRanges(mode = 'TE', results_df = results_df_local_TE)

GRanges_TE_start = GRanges_TE
end(GRanges_TE_start) = start(GRanges_TE_start)

## Genes

GRanges_gene = make_GRanges(mode = 'gene', results_df = results_df_local_gene)

GRanges_gene_extended = GRanges_gene
end(GRanges_gene_extended) = end(GRanges_gene) + 1000
start(GRanges_gene_extended) = start(GRanges_gene) - 1000

GRanges_gene_upstream = GRanges_gene
end(GRanges_gene_upstream) = start(GRanges_gene)
start(GRanges_gene_upstream) = start(GRanges_gene) - 5000

## EREs

GRanges_ERE = make_GRanges(mode = 'TE', results_df = results_df_local_ERE)

# Genes of interest

ensembl = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

AIRE_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/human_aire_dep_genes_san.csv')$ensembl_gene_id
#FEZF2_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/human_fezf2_dep_genes.csv')$ensembl_gene_id
#TRA_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/human_tra_genes.csv')$ensembl_gene_id

housekeeping_genes = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/housekeeping_genes.csv', sep = ' ')$refseq_mrna
housekeeping_genes = getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", values = housekeeping_genes, mart= ensembl)$ensembl_gene_id

#################################################################
# Permutation test
#################################################################

calculate_permutation = function(gene_sets,
                                 iterations = 100){
  
  ## 'gene sets' is a list in which each entry is a list of ENSEMBL gene IDs.

  subset = names(gene_sets)
  p_value = vector()
  Z_score = vector()

  for (i in 1:length(gene_sets)){
    
    gene_set = gene_sets[[i]]
    
    pt = permTest(A = subset(GRanges_TE_start, significant == T & log2FoldChange > 0), 
                  B = subset(GRanges_gene_extended, Geneid %in% gene_set), 
                  ntimes = iterations,
                  randomize.function = resampleRegions,
                  universe = GRanges_TE_start,
                  evaluate.function = numOverlaps,
                  alternative = 'greater',
                  verbose = TRUE)
    
    p_value[i] = pt$numOverlaps[[1]]
    Z_score[i] = pt$numOverlaps[[6]]
   
    print(pt$numOverlaps[[6]])
        
  }
  
  p_value = p.adjust(p_value, method = 'bonferroni')
  
  output = data.frame(subset = subset, 
                      p_value = p_value, 
                      Z_score = Z_score)
  
  return(output)
  
}

output = calculate_permutation(gene_sets = gene_sets)

pt = permTest(A = subset(GRanges_TE, significant == T & log2FoldChange > 0), 
              B = subset(GRanges_gene, Geneid %in% housekeeping_genes), 
              ntimes = 100,
              randomize.function = resampleRegions,
              universe = GRanges_TE,
              evaluate.function = meanDistance,
              alternative = 'less',
              verbose = TRUE)



## 

perm_test_output_A = run_perm_test(group_A = list(AIRE = subset(GRanges_gene_extended, Geneid %in% AIRE_genes),
                                             housekeeping = subset(GRanges_gene_extended, Geneid %in% housekeeping_genes),
                                             TRA = subset(GRanges_gene_extended, Geneid %in% TRA_genes),
                                             FEZF2 = subset(GRanges_gene_extended, Geneid %in% FEZF2_genes)),
                                   group_B = list(TE_up = subset(GRanges_TE_start, significant == T & log2FoldChange > 0), 
                                                  TE_down = subset(GRanges_TE_start, significant == T & log2FoldChange < 0)),
                                   universe = GRanges_gene_extended,
                                   mode = 'overlap',
                                   iterations = 1000)

perm_test_output_C = run_perm_test(group_A = list(gene_up = subset(GRanges_gene_extended, significant == T & log2FoldChange > 0), 
                                                  gene_down = subset(GRanges_gene_extended, significant == T & log2FoldChange < 0)),
                                   group_B = list(TE_up = subset(GRanges_TE_start, significant == T & log2FoldChange > 0), 
                                                  TE_down = subset(GRanges_TE_start, significant == T & log2FoldChange < 0)),
                                   universe = GRanges_gene_extended,
                                   mode = 'overlap',
                                   iterations = 1000)

perm_test_output_D = run_perm_test(group_A = list(gene_up = subset(GRanges_gene_extended, significant == T & log2FoldChange > 0), 
                                                  gene_down = subset(GRanges_gene_extended, significant == T & log2FoldChange < 0)),
                                   group_B = list(LTR = subset(GRanges_TE_start, class == 'LTR'), 
                                                  SINE = subset(GRanges_TE_start, class == 'SINE'),
                                                  LINE = subset(GRanges_TE_start, class == 'LINE'),
                                                  DNA = subset(GRanges_TE_start, class == 'DNA')),
                                   universe = GRanges_gene_extended,
                                   mode = 'overlap',
                                   iterations = 100)



## Heatmap

input = perm_test_output_A

labels = input$p_value

for (i in 1:length(labels)){
  
  old_value = labels[i]
  p = input[[2]][i]

  significance_threshold = 0.05
  
  if (p < significance_threshold){
    
    p_value = as.character(signif(p, digits = 2))
    annotation = glue('p = {p_value}')
    
  }
  
  else{
    
    annotation = 'ns'
    
  }
  
  new_value = annotation
  labels[i] = new_value
  
}

my_heatmap = pheatmap(mat = input[[1]],
                      color = colorRampPalette(rev(c("#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB")))(100),
                      cluster_rows=FALSE,
                      show_rownames=TRUE, 
                      cluster_cols=FALSE,
                      fontsize_number = 15,
                      fontsize = 15,
                      display_numbers = labels,
                      number_color = 'black',
                      border_color = 'black',
                      angle_col = '0'),
                      labels_row = c('  Up', '   - ', '  Down'),
                      labels_col = c('Up', '-', 'Down'))

save_pheatmap_png(my_heatmap, "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/geneofinterest_TE_overlap.png")

## Overlap between TSSs and TEs

df_all = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_all.csv') 
GRanges_TSS = makeGRangesFromDataFrame(df_all, keep.extra.columns = T)

df_high = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_high.csv') 
GRanges_TSS_high = makeGRangesFromDataFrame(df_high, keep.extra.columns = T)
start(GRanges_TSS_high_extended) = GenomicRanges::start(GRanges_TSS_high) - 1000
end(GRanges_TSS_high_extended) = GenomicRanges::start(GRanges_TSS_high) + 1000

df_low = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_low.csv') 
GRanges_TSS_low = makeGRangesFromDataFrame(df_low, keep.extra.columns = T)

perm_test_output_B = run_perm_test(group_A = list(TE_up = GRanges_TE_up_extended, TE_unchanged = GRanges_TE_unchanged_extended, TE_down = GRanges_TE_down),
                                   group_B = list(tss_high = GRanges_TSS_high, tss_low = GRanges_TSS_low),
                                   universe = GRanges_TE,
                                   mode = 'overlap',
                                   iterations = 1000)

my_heatmap = pheatmap(mat = perm_test_output_A[[1]],
                      color = colorRampPalette(rev(c("#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB")))(100),
                      cluster_rows=FALSE,
                      show_rownames=TRUE, 
                      cluster_cols=FALSE,
                      display_numbers = labels,
                      fontsize_number = 15,
                      fontsize = 15,
                      number_color = 'black',
                      border_color = 'black',
                      angle_col = '0',
                      labels_row = c('  Up', '   - ', '  Down'),
                      labels_col = c('Up', '-', 'Down'))


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


gene_sets = list('up_genes' = subset(GRanges_gene, significant == T & log2FoldChange > 0)$Geneid,
                 'down_genes' = subset(GRanges_gene, significant == T & log2FoldChange < 0)$Geneid)

output = calculate_odds_ratio_of_overlap(gene_sets = gene_sets,
                                         TE_set = 'up',
                                         TE_coordinates = GRanges_TE_start,
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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/odds_ratio_overlap.png", 
       width = 5, height = 5, units = "in")

up_genes = subset(GRanges_gene, significant == T & log2FoldChange > 0)$Geneid

length(TRA_genes %in% up_genes) / length(up_genes) * 100
length(housekeeping_genes %in% up_genes) / length(up_genes) * 100


#gene_sets = list('up_genes' = subset(GRanges_gene_extended, significant == T & log2FoldChange > 0),
#                 'unchanged_genes' = subset(GRanges_gene_extended, significant == F),
#                 'down_genes' = subset(GRanges_gene_extended, significant == T & log2FoldChange < 0))

gene_sets = list('AIRE_genes' = subset(GRanges_gene_extended, Geneid %in% AIRE_genes),
                 'FEZF2_genes' = subset(GRanges_gene_extended, Geneid %in% FEZF2_genes),
                 'other_genes' = subset(GRanges_gene_extended, !(Geneid %in% AIRE_genes) & !(Geneid %in% FEZF2_genes)))


TE_sets = list('up_TEs' = subset(GRanges_TE_start, significant == T & log2FoldChange > 0),
               'unchanged_TEs' = subset(GRanges_TE_start, significant == F),
               'down_TEs' = subset(GRanges_TE_start, significant == T & log2FoldChange < 0))

for (a in 1:length(TE_sets)){
  
  total_overlaps = length(findOverlaps(query = TE_sets[[a]],
                                       subject = GRanges_gene_extended))
  
  percent = vector()
  for (b in 1:length(gene_sets)){
    
    subset_overlap = length(findOverlaps(query = TE_sets[[a]],
                                         subject = gene_sets[[b]]))
    
    percent[b] = subset_overlap / total_overlaps * 100
    
  }
  
  output = data.frame(gene_set = names(gene_sets),
                      percent_overlap = percent,
                      TE_set = names(TE_sets)[a])
  
  if (a == 1){
    
    final_output = output
    
  }
  
  else{
    
    final_output = bind_rows(final_output, output)
    
  }
  
}

final_output$gene_set = factor(final_output$gene_set, levels = c('down_genes', 'unchanged_genes', 'up_genes'))
final_output$TE_set = factor(final_output$TE_set, levels = c('up_TEs', 'unchanged_TEs', 'down_TEs'))

bar_chart = ggplot(final_output, aes(x = TE_set, y = percent_overlap, fill = gene_set)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  xlab('') +
  ylab('Fraction of overlap events') +
  labs(fill= "") +
  scale_x_discrete(labels = c('Up TEs', 'Unchanged TEs', 'Down TEs')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  guides(fill = guide_legend(reverse = TRUE))

bar_chart + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                               plot.subtitle = element_text(size = 14),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text.x = element_text(size = 13, margin = margin(t = 6)),
                               axis.text.y = element_text(size = 14),
                               axis.title.y = element_text(size = 14),
                               axis.title.x = element_text(size = 14, margin = margin(t = 6)),
                               axis.line = element_line(size = 0.8),
                               panel.border = element_blank(),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 14),
                               legend.position="top")

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/fraction_of_overlaps_with_GOIs.png", 
       width = 5, height = 5, units = "in")

#################################################################
# ereMAPs
#################################################################

ereMAPs = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/ereMAPs_Larouche.csv', header = T)
ereMAPs$Start <- as.numeric(gsub(",","",ereMAPs$Start))
ereMAPs$End <- as.numeric(gsub(",","",ereMAPs$End))

GRanges_ereMAPs = makeGRangesFromDataFrame(ereMAPs, keep.extra.columns = T)

TE_ereMAP_overlaps = findOverlaps(query = GRanges_ERE,
                                 subject = GRanges_ereMAPs)

TE_ereMAP_overlaps = data.frame(TE_ereMAP_overlaps)

validated = vector()
ereMAP_loci = vector()
for (entry in 1:nrow(TE_ereMAP_overlaps)){
  
  TE_hit = TE_ereMAP_overlaps[entry, ]$queryHits
  ereMAP_hit = TE_ereMAP_overlaps[entry, ]$subjectHits
  
  if (GRanges_ERE[TE_hit, ]$gene == GRanges_ereMAPs[ereMAP_hit, ]$ERE.family){
    
    validated[entry] = T
    
  }
  
  else{
    
    validated[entry] = F

  }
  
  ereMAP_loci[entry] = GRanges_ERE[TE_hit, ]$ID
  
}

results_df_local_ERE = mutate(results_df_local_ERE, ereMAP = case_when(ID %in% ereMAP_loci ~ T,
                                                                      !(ID %in% ereMAP_loci) ~ F))

output = generate_contingency(input = results_df_local_ERE,
                              condition_A = list('up_TEs' = 'significant == T & log2FoldChange > 0',
                                                 'down_TEs' = 'significant == T & log2FoldChange < 0',
                                                 'unchanged_TEs' = 'significant == F'),
                              condition_B = list('ereMAP' = 'ereMAP == T', 'not_ereMAP' = 'ereMAP == F'))

vcd::mosaic(~condition_A+condition_B, data = output, direction = c('v', 'h'), shade = T)

#################################################################
# Ranked by expression
#################################################################

normalized_counts = data.frame(assay(vs_dds_local))

normalized_counts_ERE = extract_subset(mode = 'ERE', input = normalized_counts)

normalized_counts_ERE_hi = dplyr::select(normalized_counts_ERE, 
                                     c('pt226_hi_fastp_1.fastq_Aligned.out.bam', 
                                       'pt221_lo_fastp_1.fastq_Aligned.out.bam', 
                                       'pt214_lo_fastp_1.fastq_Aligned.out.bam'))

rowMeans(normalized_counts_hi)
