#################################################################
# Libraries
#################################################################

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

library(GenomicRanges)
library(regioneR)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(glue)

#################################################################
# Functions
#################################################################

build_count_table = function(dds, results_df, group, mode, by){
  
  normalized_counts = as.data.frame(counts(dds, normalized = TRUE))
  
  normalized_counts = cbind(ID = rownames(normalized_counts), normalized_counts)
  
  ## Merge with results_df, determine mean and remove unnecessary columns
  
  normalized_counts = merge(normalized_counts, results_df, by = 'ID')
  normalized_counts = select(normalized_counts, -c('pt214_mTEC-lo_new', 'pt221_mTEC-lo_new', 'pt226_mTEC-lo_new'))
  normalized_counts$mean = rowMeans(normalized_counts[2:4])
  normalized_counts = select(normalized_counts, -c('pt214_mTEC-hi_new', 'pt221_mTEC-hi_new', 'pt226_mTEC-hi_new'))
  
  normalized_counts = mutate(normalized_counts, class = sub("\\?", "", class))
  
  ## Input is determined by subsetting normalized_counts by logical conditions specified in 'group'
  
  for (i in group){
    
    if (i == 'all'){
      
      input = normalized_counts
      
    }
    
    if (i == 'diff_regulated'){
      
      sigGenes = results_df[results_df$significant == TRUE,]$ID
      input = normalized_counts[normalized_counts$ID %in% sigGenes,]
      
    }
    
    if (i == 'up_regulated'){
      
      upGenes = results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange > 0), ]$ID
      input = normalized_counts[normalized_counts$ID %in% upGenes,]
      
    }
    
    if (i == 'down_regulated'){
      
      downGenes = results_df[(results_df$significant == TRUE) & (results_df$log2FoldChange < 0), ]$ID
      input = normalized_counts[normalized_counts$ID %in% downGenes,]
      
    }
    
    if (i == 'diff_regulated_both'){
      
      subset = results_df[results_df$significant == TRUE,]$ID
      input = normalized_counts[rownames(normalized_counts) %in% subset,]
      
      subset_2 = results_df_transcripts[results_df_transcripts$significant == TRUE,]$gene
      input = filter(input, gene %in% subset_2)
      
    }
    
    ########
    
    if (mode == 'class'){
      
      if (by == 'normalized_reads'){
        
        input = group_by(input, class) %>% summarize(summary = sum(mean))
        
      }
      
      if (by == 'locus'){
        
        input = group_by(input, class) %>% summarize(summary = n())
        
      }
      
    }
    
    if (mode == 'LTR_family'){
      
      if (by == 'normalized_reads'){
        
        input = input %>% 
          filter(class == 'LTR') %>%
          group_by(family) %>% 
          summarize(summary = sum(mean))
        
      }
      
      if (by == 'locus'){
        
        input = input %>% 
          filter(class == 'LTR') %>%
          group_by(family) %>% 
          summarize(summary = n())
        
      }
      
    }
    
    if (mode == 'overlap'){
      
      if (by == 'normalized_reads'){
        
        input = group_by(input, overlap_status) %>% summarize(summary = sum(mean))
        
      }
      
      if (by == 'locus'){
        
        input = group_by(input, overlap_status) %>% summarize(summary = n())
        
      }
      
    }
    
    ######
    
    input = as.data.frame(input)
    input$group = i
    
    if (match(i, group) == 1){
      
      output = input
      
    }
    
    if (match(i, group) != 1) {
      
      output = bind_rows(output, input)
      
    }  
    
  }
  
  output = output %>%
    group_by(group) %>%
    mutate(percent = summary / sum(summary) * 100)
  
  return(output)
  
}

make_GRanges = function(mode, results_df){
  
  if (mode == 'TE'){
    
    annotation = read.table(file = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt", header = 1)
    annotation = separate(annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':')
    annotation = separate(annotation, start.stop, into = c('start', 'end'), sep = '-')
    annotation = rename(annotation, locus = TE)
    
    df = merge(results_df, annotation, by = 'locus')
    
  }
  
  if (mode == 'gene'){
    
    annotation = read.table(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/annotation_tables/gencode.v38_gene_annotation_table.txt', header = 1)
    annotation = select(annotation, c('Geneid', 'Chromosome', 'Start', 'End', 'Strand', 'Class'))
    annotation = rename(annotation, chr = Chromosome, start = Start, end = End, strand = Strand)
    
    df =  merge(results_df, annotation, by = 'Geneid')
    
  }
  
  if (mode == 'TSS'){
    
    df = select(results_df, c('Chr', 'Start', 'Stop', 'Strand'))
    
  }
  
  output = makeGRangesFromDataFrame(df, keep.extra.columns = T)
  
  return(output)
  
}

run_perm_test = function(TSS_groups, TE_groups, mode = 'distance'){
  
  Z_score_output = matrix(, nrow = length(TSS_groups), ncol = length(TE_groups))
  rownames(Z_score_output) = c('High TSS', 'Low TSS')
  colnames(Z_score_output) = c('UP TE', 'Unchanged TE', 'DOWN TE')
  
  p_value_output = Z_score_output
  
  row_number = 1
  
  for (TSS_group in TSS_groups){
    
    column_number = 1
    
    for (TE_group in TE_groups){
      
      print(glue('Starting row {row_number}, column {column_number}'))
      
      if (mode == 'distance'){
        
        pt = permTest(A = TSS_group, 
                      B = TE_group, 
                      ntimes = 1000,
                      randomize.function = resampleRegions,
                      universe = GRanges_TSS_all,
                      evaluate.function = meanDistance,
                      alternative = 'less',
                      verbose = TRUE)
        
        print(pt$meanDistance[[1]])
        
        p_value_output[row_number, column_number] = pt$meanDistance[[1]]
        Z_score_output[row_number, column_number] = pt$meanDistance[[6]]
        
      }
      
      if (mode == 'overlap'){
        
        pt = permTest(A = TE_group, 
                      B = gene_group, 
                      ntimes = 2000,
                      randomize.function = resampleRegions,
                      universe = GRanges_TE,
                      evaluate.function = numOverlaps,
                      alternative = 'greater',
                      verbose = TRUE)
        
        p_value_output[row_number, column_number] = pt$numOverlaps[[1]]
        Z_score_output[row_number, column_number] = pt$numOverlaps[[6]]
        
        print(Z_score_output)
        
      }
      
      column_number = column_number + 1
      
    }
    
    row_number = row_number + 1
    
  }
  
  output = list('Z score' = Z_score_output, 'p value' = p_value_output)
  
  return(output)
  
}


#################################################################
# PCA w/ GTEx data (A)
#################################################################

pcaData = plotPCA(vs_dds_transcripts_TE, intgroup='tissue', returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

PCA = ggplot(pcaData, aes(PC1, PC2, fill = tissue)) + 
  geom_point(size=4, shape = 21, stroke = 0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  labs(fill= "Tissue") +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff')) 

PCA + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                         plot.subtitle = element_text(size = 14),
                         axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),
                         axis.title = element_text(size = 14),
                         axis.line = element_line(size = 0.8),
                         panel.border = element_blank(),
                         legend.text = element_text(size = 15),
                         legend.title = element_text(size = 0),
                         legend.position="top",
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/A_pca.png", 
       width = 5.25, height = 6, units = "in")

#################################################################
# Volcano plot (B)
#################################################################

input = results_df_transcripts_TE

input = mutate(input, ID = sub("\\?", "", ID))
input = mutate(input, class = sub("\\?", "", class))

## Version 1

input_v1 = mutate(input, grouped_class = case_when(class == 'DNA' ~ 'Other', 
                                                class == 'LINE' ~ 'Other',
                                                class == 'LTR' ~ 'LTR',
                                                class == 'SINE' ~ 'Other',
                                                class == 'Retroposon' ~ 'Other',
                                                class == 'Satellite' ~ 'Other'))

volcano_plot = ggplot() +
  geom_point(data = input_v1, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input_v1, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = grouped_class), size = 2.5, alpha = 1, shape = 21, stroke = 0) +
  geom_point(data = subset(input_v1, significant == FALSE), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 1, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(FDR)')) +
  xlim(-2, 3) +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff')) +
  labs(fill= "")

## Version 2

input_v2 = mutate(input, grouped_class = case_when(class == 'LINE' ~ 'LINE',
                                                   class == 'LTR' ~ 'LTR',
                                                   class == 'SINE' ~ 'SINE',
                                                   class == 'Retroposon' ~ 'Other',
                                                   class == 'Satellite' ~ 'Other',
                                                   class == 'RC' ~ 'Other',
                                                   class == 'DNA' ~ 'DNA',
                                                   class == 'RNA' ~ 'Other',
                                                   class == 'Unknown' ~ 'Other'))

input_v2$grouped_class = factor(input_v2$grouped_class, levels = c('LTR', 'DNA', 'LINE', 'SINE', 'Other'))

volcano_plot = ggplot() +
  geom_point(data = input_v2, aes(x = log2FoldChange, y = -log10(padj)), color = alpha('#9B9A99', 0.6)) +
  geom_point(data = subset(input_v2, significant == TRUE), aes(x = log2FoldChange, y = -log10(padj), fill = grouped_class), size = 2.5, alpha = 1, shape = 21, stroke = 0) +
  geom_point(data = subset(input_v2, significant == FALSE), aes(x = log2FoldChange, y = -log10(padj)), size = 1, alpha = 1, shape = 21, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  xlab(expression('log'[2]*'(fold-change)')) +
  ylab(expression('-log'[10]*'(adjusted p-value)')) +
  xlim(-2, 3) +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff', '#55A257', '#E93C00', '#9A9A9A')) +
  labs(fill= "")

#ggtitle('mTEC-hi vs mTEC-lo', 'TE transcripts') 

## Plot

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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/B_volcano_v2.png", 
       width = 5.25, height = 5.25, units = "in")

#################################################################
# Stacked bars (C)
#################################################################

count_table = build_count_table(dds_transcripts_TE, 
                                results_df_transcripts_TE, 
                                group = c('all', 'down_regulated', 'up_regulated'),
                                mode = 'class',
                                by = 'normalized_reads')

count_table = mutate(count_table, grouped_class = case_when(class == 'LINE' ~ 'LINE',
                                                   class == 'LTR' ~ 'LTR',
                                                   class == 'SINE' ~ 'SINE',
                                                   class == 'Retroposon' ~ 'Other',
                                                   class == 'Satellite' ~ 'Other',
                                                   class == 'RC' ~ 'Other',
                                                   class == 'DNA' ~ 'DNA',
                                                   class == 'RNA' ~ 'Other',
                                                   class == 'Unknown' ~ 'Other'))


count_table = group_by(count_table, group, grouped_class) %>% 
  summarize(percent = sum(percent))

count_table$grouped_class = factor(count_table$grouped_class, levels = c('Other', 'SINE', 'LINE', 'DNA', 'LTR'))

bar_chart = ggplot(count_table, aes(x = group, y = percent, fill = grouped_class)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  xlab('') +
  ylab('Fraction of normalized reads') +
  labs(fill= "") +
  scale_x_discrete(labels = c('All', 'Downregulated', 'Upregulated')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c('#9A9A9A', '#E93C00', '#55A257', '#dd8452ff', '#4c72b0ff'))

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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/C_stacked-bars.png", 
       width = 5, height = 5, units = "in")

## Chi-square test

contingency = pivot_wider(count_table, id_cols = grouped_class, names_from = group, values_from = sum) %>% 
  mutate_if(is.numeric, list(~replace_na(., 0))) %>% 
  as.data.frame()

rownames(contingency) = contingency$grouped_class
contingency = dplyr::select(contingency, -grouped_class)

chisq.test(contingency)

chisq.posthoc.test::chisq.posthoc.test(contingency)

#################################################################
# Genomic position analysis (D)
#################################################################

GRanges_TE = make_GRanges(mode = 'TE',
                          results_df = results_df_local_TE)
GRanges_TE_start = GRanges_TE
end(GRanges_TE_start) = GenomicRanges::start(GRanges_TE_start) + 100

GRanges_TE_up = make_GRanges(mode = 'TE',
                             results_df = results_df_local_TE_up)
GRanges_TE_up_start = GRanges_TE_up
end(GRanges_TE_up_start) = GenomicRanges::start(GRanges_TE_up_start) + 100

GRanges_TE_down = make_GRanges(mode = 'TE',
                             results_df = results_df_local_TE_down)
GRanges_TE_down_start = GRanges_TE_down
end(GRanges_TE_down_start) = GenomicRanges::start(GRanges_TE_down_start) + 100

GRanges_TE_unchanged = make_GRanges(mode = 'TE',
                             results_df = results_df_local_TE_unchanged)
GRanges_TE_unchanged_start = GRanges_TE_unchanged
end(GRanges_TE_unchanged_start) = GenomicRanges::start(GRanges_TE_unchanged_start) + 100

df_high = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_high.csv')
GRanges_TSS_high = make_GRanges(mode = 'TSS',
                                results_df = df_high)

df_low = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_low.csv')
GRanges_TSS_low = make_GRanges(mode = 'TSS',
                                results_df = df_low)


df_all = read.csv(file = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/TSS/df_all.csv')
GRanges_TSS_all = make_GRanges(mode = 'TSS',
                               results_df = df_all)

output = run_perm_test(TSS_groups = list(GRanges_TSS_high, GRanges_TSS_low),
                       TE_groups = list(GRanges_TE_up_start, GRanges_TE_unchanged_start, GRanges_TE_down_start))

#################################################################
# Fraction of reads mapping to TEs (supplement A?)
#################################################################

#################################################################
#  (supplement B?)
#################################################################

count_table = build_count_table(dds_transcripts_TE, 
                                results_df_transcripts_TE, 
                                group = c('all', 'down_regulated', 'up_regulated'),
                                mode = 'LTR_family',
                                by = 'normalized_reads')

count_table = mutate(count_table, grouped_class = case_when(family == 'ERV1' ~ 'ERV1',
                                                            family == 'ERVK' ~ 'ERVK',
                                                            family == 'ERVL' ~ 'ERVL',
                                                            family == 'ERVL-MaLR' ~ 'ERVL-MaLR',
                                                            family == 'Gypsy' ~ 'Other',
                                                            family == 'LTR' ~ 'Other'))

count_table = group_by(count_table, group, grouped_class) %>% 
  summarize(percent = sum(percent))

count_table$grouped_class = factor(count_table$grouped_class, levels = c('Other', 'ERVL-MaLR', 'ERVL', 'ERVK', 'ERV1'))

bar_chart = ggplot(count_table, aes(x = group, y = percent, fill = grouped_class)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  xlab('') +
  ylab('Fraction of normalized reads') +
  labs(fill= "") +
  scale_x_discrete(labels = c('All', 'Downregulated', 'Upregulated')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c('#9A9A9A', '#E93C00', '#55A257', '#dd8452ff', '#4c72b0ff'))

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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/SB_stacked-bars.png", 
       width = 6, height = 5, units = "in")