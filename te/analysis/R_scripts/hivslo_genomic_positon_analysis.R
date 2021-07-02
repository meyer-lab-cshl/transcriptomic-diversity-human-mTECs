library(GenomicRanges)
library(regioneR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(glue)

#################################################################
# Functions
#################################################################

make_GRanges = function(mode, results_df){
  
  if (mode == 'TE'){
    
    annotation = read.table(file = "/Users/mpeacey/TE_thymus/analysis/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt", header = 1)
    annotation = separate(annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':')
    annotation = separate(annotation, start.stop, into = c('start', 'end'), sep = '-')
    annotation = rename(annotation, locus = TE)
    
    df = merge(results_df, annotation, by = 'locus')
    
  }
  
  if (mode == 'gene'){
    
    annotation = read.table(file = '/Users/mpeacey/TE_thymus/analysis/annotation_tables/gencode.v38_gene_annotation_table.txt', header = 1)
    annotation = select(annotation, c('Geneid', 'Chromosome', 'Start', 'End', 'Strand', 'Class'))
    annotation = rename(annotation, chr = Chromosome, start = Start, end = End, strand = Strand)
    
    df =  merge(results_df, annotation, by = 'Geneid')
    
  }
  
  output = makeGRangesFromDataFrame(df, keep.extra.columns = T)
  
  return(output)
  
}

run_perm_test = function(gene_groups, TE_groups, mode = 'overlap'){
  
  Z_score_output = matrix(, nrow = length(gene_groups), ncol = length(TE_groups))
  rownames(Z_score_output) = c('UP gene', 'Unchanged gene', 'DOWN gene')
  colnames(Z_score_output) = c('UP TE', 'Unchanged TE', 'DOWN TE')
  
  p_value_output = Z_score_output

  row_number = 1
  
  for (gene_group in gene_groups){
    
    column_number = 1
    
    for (TE_group in TE_groups){
      
      print(glue('Starting row {row_number}, column {column_number}'))
      
      if (mode == 'distance'){
        
        pt = permTest(A = TE_group, 
                      B = gene_group, 
                      ntimes = 1000,
                      randomize.function = resampleRegions,
                      universe = GRanges_TE,
                      evaluate.function = meanDistance,
                      alternative = 'less',
                      verbose = TRUE)
        
        p_value = pt$meanDistance[[1]]
        Z_score = pt$meanDistance[[6]]
        
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

correlate_fold_change = function(query, subject){
  
  query_fold_change = rep(NA, length(query))
  subject_fold_change = rep(NA, length(query))
  
  for (i in 1:length(query)){
    print(i)
    query_fold_change[i] = query[i]$log2FoldChange[1]
    
    overlapping_subjects = findOverlaps(query = query[i],
                                        subject = subject)
    
    overlapping_subjects = as.data.frame(overlapping_subjects)
    hit_indices = overlapping_subjects$subjectHits
    
    fold_changes = rep(NA, length(hit_indices))
    index_number = 1
    
    for (index in hit_indices){
      
      fold_changes[index_number] = 2 ** subject[index]$log2FoldChange[1]
      index_number = index_number + 1
      
    }
    
    subject_fold_change[i] = log2(mean(fold_changes))
    
  }
  
  output = data.frame('query_fold_change' = query_fold_change, 'subject_fold_change' = subject_fold_change)
  
  return(output)
  
}

build_count_table = function(group){
  
  ## Input is determined by subsetting normalized_counts by logical conditions specified in 'group'
  
  for (i in group){
    
    if (i == 'all'){
      
      input = as.data.frame(GRanges_TE_annotated)
      
    }
    
    if (i == 'diff_regulated'){
      
      input = as.data.frame(GRanges_TE_annotated) %>%
        dplyr::filter(significant == TRUE)
      
    }
    
    if (i == 'up_regulated'){
      
      input = as.data.frame(GRanges_TE_annotated) %>%
        dplyr::filter((significant == TRUE) & (log2FoldChange > 0))
      
    }
    
    input = group_by(input, Class) %>% summarize(summary = n())
    
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

#################################################################
# GRanges
################################################################# 

## Make GRanges_gene objects

GRanges_gene = make_GRanges(mode = 'gene',
                            results_df = results_df_local_gene)
saveRDS(GRanges_gene, "~/TE_thymus/analysis/cluster/objects/GRanges_gene.rds")

GRanges_gene_up = make_GRanges(mode = 'gene',
                               results_df = results_df_local_gene_up)
saveRDS(GRanges_gene_up, "~/TE_thymus/analysis/cluster/objects/GRanges_gene_up.rds")

GRanges_gene_unchanged = make_GRanges(mode = 'gene',
                               results_df = results_df_local_gene_unchanged)
GRanges_gene_down = make_GRanges(mode = 'gene',
                               results_df = results_df_local_gene_down)
saveRDS(GRanges_gene_down, "~/TE_thymus/analysis/cluster/objects/GRanges_gene_down.rds")

GRanges_gene_sigdiff = make_GRanges(mode = 'gene',
                                    results_df = results_df_local_gene_sigdiff)
saveRDS(GRanges_gene_sigdiff, "~/TE_thymus/analysis/cluster/objects/GRanges_gene_sigdiff.rds")

## Make GRanges_TE objects

GRanges_TE = make_GRanges(mode = 'TE',
                          results_df = results_df_local_TE)
saveRDS(GRanges_TE, "~/TE_thymus/analysis/cluster/objects/GRanges_TE.rds")

GRanges_TE_up = make_GRanges(mode = 'TE',
                               results_df = results_df_local_TE_up)
saveRDS(GRanges_TE, "~/TE_thymus/analysis/cluster/objects/GRanges_TE_up.rds")

GRanges_TE_unchanged = make_GRanges(mode = 'TE',
                                      results_df = results_df_local_TE_unchanged)
saveRDS(GRanges_TE, "~/TE_thymus/analysis/cluster/objects/GRanges_TE_unchanged.rds")

GRanges_TE_down = make_GRanges(mode = 'TE',
                                 results_df = results_df_local_TE_down)
saveRDS(GRanges_TE, "~/TE_thymus/analysis/cluster/objects/GRanges_TE_down.rds")

GRanges_TE_sigdiff = make_GRanges(mode = 'TE',
                                    results_df = results_df_local_TE_sigdiff)
saveRDS(GRanges_TE_sigdiff, "~/TE_thymus/analysis/cluster/objects/GRanges_TE_sigdiff.rds")

## Annotate by overlap

annotate_features = function(input, mode){
  
  if (mode == 'hannah'){
    
    features = readRDS('~/TE_thymus/analysis/features.rds')
    
    annotated_GRanges_TE = sort(input)
    
    annotated_GRanges_TE$TSS_10_genes <- overlapsAny(input, features$TSS10_genes)
    annotated_GRanges_TE$TSS_10_trans <- overlapsAny(input, features$TSS10_transcripts)
    annotated_GRanges_TE$TSS_100_genes <- overlapsAny(input, features$TSS100_genes)
    annotated_GRanges_TE$TSS_100_trans <- overlapsAny(input, features$TSS100_transcripts)
    annotated_GRanges_TE$first_exon <- overlapsAny(input, features$first_exon)
    annotated_GRanges_TE$other_exon <- overlapsAny(input, features$other_exon)
    annotated_GRanges_TE$intron <- overlapsAny(input, features$introns)
    annotated_GRanges_TE$TTS_100_genes <- overlapsAny(input, features$TTS100_genes)
    annotated_GRanges_TE$TTS_200_genes <- overlapsAny(input, features$TTS200_genes)
    annotated_GRanges_TE$TTS_100_trans <- overlapsAny(input, features$TTS100_trans)
    annotated_GRanges_TE$TTS_200_trans <- overlapsAny(input, features$TTS200_trans)
    annotated_GRanges_TE$downstream_gene <- overlapsAny(input, features$downstream_genes)
    annotated_GRanges_TE$upstream_gene <- overlapsAny(input, features$upstream_genes)
    annotated_GRanges_TE$antisense <- overlapsAny(input, features$antisense)
    
  }
  
  if (mode == 'matty'){
    
    annotation = read.table(file = '/Users/mpeacey/TE_thymus/analysis/annotation_tables/gencode.v38_gene_annotation_table.txt', header = 1)
    annotation = rename(annotation, chr = Chromosome, start = Start, end = End, strand = Strand)
    
    new_annotation = merge(annotation, results_df_local_gene, by = 'Geneid', all = T)
    
    features = makeGRangesFromDataFrame(new_annotation, keep.extra.columns = T)

    annotated_GRanges_TE = splicejam::annotateGRfromGR(GR1 = sort(input), GR2 = features)
      
  
  }
  
  return(annotated_GRanges_TE)
  
}

GRanges_TE_annotated = annotate_features(input = GRanges_TE, mode = 'matty')

count_table = build_count_table(group = c('all', 'diff_regulated', 'up_regulated')) %>% dplyr::filter(percent > 1)

## Plot stacked bars

bar_chart = ggplot(count_table, aes(x = group, y = percent, fill = Class)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  scale_fill_brewer(palette = "Set1") +
  xlab('') +
  ylab('Fraction of normalized reads') +
  ggtitle('TEs in mTEC-HI cells', 'Subset by physical overlap with genes')

bar_chart + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                               plot.subtitle = element_text(size = 14),
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_line(color = 'grey'),
                               panel.grid.minor.y = element_blank(),
                               axis.text.x = element_text(size = 13),
                               axis.text.y = element_text(size = 14),
                               axis.title = element_text(size = 14),
                               axis.line = element_line(size = 0.8),
                               panel.border = element_blank(),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 14))

#ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/21-03-25/slide_3.png", 
#       width = 20, height = 13, units = "cm")

#################################################################
# regioneR
#################################################################

GRanges_gene_up_extended = GRanges_gene_up
start(GRanges_gene_up_extended) = GenomicRanges::start(GRanges_gene_up) - 2000
end(GRanges_gene_up_extended) = GenomicRanges::end(GRanges_gene_up) + 2000

GRanges_gene_down_extended = GRanges_gene_down
start(GRanges_gene_down_extended) = GenomicRanges::start(GRanges_gene_down) - 2000
end(GRanges_gene_down_extended) = GenomicRanges::end(GRanges_gene_down) + 2000

GRanges_gene_unchanged_extended = GRanges_gene_unchanged
start(GRanges_gene_unchanged_extended) = GenomicRanges::start(GRanges_gene_unchanged) - 2000
end(GRanges_gene_unchanged_extended) = GenomicRanges::end(GRanges_gene_unchanged) + 2000

gene_groups = list(GRanges_gene_up_extended, GRanges_gene_unchanged_extended, GRanges_gene_down_extended)
TE_groups = list(GRanges_TE_up, GRanges_TE_unchanged, GRanges_TE_down)

saveRDS(gene_groups, "~/TE_thymus/analysis/cluster/objects/gene_groups.rds")
saveRDS(TE_groups, "~/TE_thymus/analysis/cluster/objects/TE_groups.rds")

output = run_perm_test(gene_groups, TE_groups, mode = 'overlap')

## Heatmap

my_heatmap = pheatmap(mat = output[[1]], 
                      cluster_rows=FALSE,
                      show_rownames=TRUE, 
                      cluster_cols=FALSE,
                      display_numbers = T,
                      fontsize_number = 15,
                      number_color = 'black')

 color = colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(100),



my_heatmap = pheatmap(mat = output, 
                      cluster_rows=FALSE,
                      show_rownames=TRUE, 
                      cluster_cols=FALSE,
                      color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                "RdYlBu")))(30),
                      display_numbers = T,
                      fontsize_number = 15,
                      number_color = 'black',
                      breaks = seq(-150, 150, by = 10))

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "/Users/mpeacey/TE_thymus/analysis/Plots/TE_local/hi_vs_lo_genomic_poisitons_heatmap_overlap.png")


#################################################################
# Correlation method
#################################################################

## query: GRanges_gene, subject: GRanges_TE

output = readRDS("~/TE_thymus/analysis/cluster/gene_vs_TE.rds")

## query: GRanges_gene_sigdiff, subject: GRanges_TE_sigdiff

output = correlate_fold_change(query = GRanges_gene_sigdiff, subject = GRanges_TE_sigdiff)

correlation = ggplot(data = output, aes(x = query_fold_change, y = subject_fold_change)) + 
  geom_point(alpha = 0.6, size = 0.5) +
  geom_smooth(method='lm') +
  xlab(expression('Log'[2]*' FC (differentially expressed genes)')) +
  ylab(expression('Log'[2]*' mean FC (overlapping differentially expressed TEs)')) +
  ggtitle('Correlated expression of overlapping genes and TEs', 'Differentially expressed genes/TEs only')

## query: GRanges_TE_sigdiff, subject: GRanges_gene_sigdiff

output = correlate_fold_change(query = GRanges_TE_sigdiff, subject = GRanges_gene_sigdiff)

correlation_test = cor.test(x = output$query_fold_change, y = output$subject_fold_change, method = 'pearson')

r_squared = as.character(correlation_test$estimate ** 2)

correlation = ggplot(data = output, aes(x = subject_fold_change, y = query_fold_change)) + 
  geom_point(alpha = 0.6, size = 0.5) +
  geom_smooth(method='lm') +
  ylab(expression('Log'[2]*' FC (differentially expressed TEs)')) +
  xlab(expression('Log'[2]*' mean FC (overlapping differentially expressed genes)')) +
  ggtitle('Correlated expression of overlapping genes and TEs', 'Differentially expressed genes/TEs only')

## Plot and test

correlation + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                               plot.subtitle = element_text(size = 14),
                               axis.text.x = element_text(size = 14),
                               axis.text.y = element_text(size = 14),
                               axis.title = element_text(size = 14),
                               axis.line = element_line(size = 0.8),
                               panel.border = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/TE_local/correlated_expression_sigdiff_only.png", 
       width = 25, height = 15, units = "cm")


###

distanceToNearest(GRanges_TE, GRanges_gene_up) %>%
  as.data.frame() %>%
  ggplot(aes(x = distance)) + 
  geom_histogram(binwidth = 1000) +
  xlim(-1000, 250000)


