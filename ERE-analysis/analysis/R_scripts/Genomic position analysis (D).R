library(DESeq2)
library(tidyverse)
library(GenomicRanges)
library(regioneR)
library(pheatmap)
library(RColorBrewer)
library(glue)

functions_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_functions/"
functions = c('extract_subset', 'differential_expression', 'process_DESeq2_results', 'make_GRanges', 'run_perm_test')

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
                        contrast = c('tissue', 'mTEC-hi', 'mTEC-lo'))

results_local_gene = extract_subset(mode = 'gene', input = results_local)
results_local_TE = extract_subset(mode = 'TE', input = results_local)

results_df_local_gene = process_DESeq2_results(results = results_local_gene, mode = 'Gene')
results_df_local_TE = process_DESeq2_results(results = results_local_TE, mode = 'TE_local')

results_df_local_gene_up = filter(results_df_local_gene, (significant == T) & (log2FoldChange > 0))
results_df_local_gene_unchanged = filter(results_df_local_gene, significant == F)
results_df_local_gene_down = filter(results_df_local_gene, (significant == T) & (log2FoldChange < 0))
results_df_local_gene_sigdiff = filter(results_df_local_gene, significant == T)

results_df_local_TE_up = filter(results_df_local_TE, (significant == T) & (log2FoldChange > 0))
results_df_local_TE_unchanged = filter(results_df_local_TE, significant == F)
results_df_local_TE_down = filter(results_df_local_TE, (significant == T) & (log2FoldChange < 0))
results_df_local_TE_sigdiff = filter(results_df_local_TE, significant == T)

#################################################################
# Genomic position analysis (D)
#################################################################

## Generate GRanges objects

GRanges_TE = make_GRanges(mode = 'TE',
                          results_df = results_df_local_TE)
GRanges_TE_start = GRanges_TE
end(GRanges_TE_start) = GenomicRanges::start(GRanges_TE_start) + 100

GRanges_TE_up = make_GRanges(mode = 'TE',
                             results_df = results_df_local_TE_up)
GRanges_TE_up_start = GRanges_TE_up
end(GRanges_TE_up_start) = GenomicRanges::start(GRanges_TE_up) + 100
start(GRanges_TE_up_start) = GenomicRanges::start(GRanges_TE_up) - 100

GRanges_TE_down = make_GRanges(mode = 'TE',
                               results_df = results_df_local_TE_down)
GRanges_TE_down_start = GRanges_TE_down
end(GRanges_TE_down_start) = GenomicRanges::start(GRanges_TE_down_start) + 100

GRanges_TE_unchanged = make_GRanges(mode = 'TE',
                                    results_df = results_df_local_TE_unchanged)
GRanges_TE_unchanged_start = GRanges_TE_unchanged
end(GRanges_TE_unchanged_start) = GenomicRanges::start(GRanges_TE_unchanged_start) + 100

## A: Overlap between genes and TEs

perm_test_output_A = run_perm_test(group_A = list(TE_up = GRanges_TE_up, TE_unchanged = GRanges_TE_unchanged, TE_down = GRanges_TE_down),
                                   group_B = list(gene_up = GRanges_gene_up, gene_unchanged = GRanges_gene_unchanged, gene_down = GRanges_gene_down),
                                   universe = GRanges_TE,
                                   mode = 'overlap',
                                   iterations = 1000)

labels = signif(perm_test_output_A$evaluation, digits = 3)

for (i in 1:length(labels)){
  
  old_value = labels[i]
  p = perm_test_output_A[[2]][i]
  print(p)
  
  significance_threshold = 0.05 / length(labels)
  
  if (p < significance_threshold){
    
    p_value = as.character(signif(p, digits = 1))
    annotation = glue('p = {p_value}')
    
  }
  
  else{
    
    annotation = 'ns'
    
  }
  
  new_value = glue('{old_value}%, {annotation}')
  labels[i] = new_value
  
}

print(labels)

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

save_pheatmap_png(my_heatmap, "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/gene_TE_overlap.png")

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
TE.counts = subset(TE.counts, class == 'LTR' | class == 'LINE')

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
prefix = "EREs"

lm = TEffectR::apply_lm(gene.annotation = gene.annotation,
                       gene.counts = gene.counts,
                       repeat.counts = sum.repeat.counts,
                       covariates = covariates,
                       prefix = prefix)

lm_results = read.table("~/Desktop/thymus-epitope-mapping/ERE-analysis/LTRs -lm-results.tsv", header= T, sep="\t") %>%
  mutate(adjusted.p.value = p.adjust(model.p.value, method = 'BH', n = nrow(lm_results))) %>%
  mutate(significant = case_when(adjusted.p.value < 0.05 ~ T,
                                 adjusted.p.value >= 0.05 ~ F))

volcano_plot = ggplot() +
  geom_point(data = lm_results, aes(x = r.squared, y = -log10(adjusted.p.value)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(lm_results, significant == TRUE), aes(x = r.squared, y = -log10(adjusted.p.value), fill = significant), size = 2.5, alpha = 1, shape = 21, stroke = 0) +
  geom_point(data = subset(lm_results, significant == FALSE), aes(x = r.squared, y = -log10(adjusted.p.value)), size = 1, alpha = 1, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('R squared')) +
  ylab(expression('-log'[10]*'(p-value)'))

# ggrepel::geom_label_repel(data = subset(lm_results, significant == TRUE), aes(x = r.squared, y = -log10(model.p.value), label = GeneName))

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
