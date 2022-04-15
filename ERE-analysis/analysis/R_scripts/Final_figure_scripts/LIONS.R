library(biomaRt)
library(GenomicRanges)
library(ggvenn)
library(tidyverse)

working_directory = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis'

#################################################################
#  Data import
#################################################################

## From Chimera 

Chimera = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/External_packages/LIONS/mTEC-analysis/mTEC-analysis.chimera', 
                   sep = '\t', header = T)
Chimera['LO_occurence'] = Chimera['Normal_Occ']
Chimera['HI_occurence'] = Chimera['Cancer_Occ']
Chimera = select(Chimera, -c('Normal_Occ', 'Cancer_Occ'))

Chimera = mutate(Chimera, total_occurence = LO_occurence + HI_occurence) %>%
  separate(col = repeatName, into = c('gene', 'class', 'family'), sep = ':', remove = F)

## Remove non EREs

Chimera = Chimera[Chimera$class != 'DNA', ]

## Remove un-needed columns

Chimera = Chimera[, c('repeatName', 'gene', 'class', 'family', 
                      'coordinates', 'RStrand', 'Contribution', 'RepeatID', 
                      'LIBRARY', 'Geneid', 'Gene_strand', 'Gene_type', 
                      'LO_occurence', 'HI_occurence', 'total_occurence')]

## Filtering Chimeric-gene transcription events

sense_Chimera = subset(Chimera, !is.na(Geneid))

for (entry in 1:nrow(sense_Chimera)){
  
  repeat_strand = sense_Chimera[entry, 'RStrand']
  
  gene_list = unlist(stringr::str_split(sense_Chimera[entry, 'Geneid'], ';'))
  gene_strand_list = unlist(stringr::str_split(sense_Chimera[entry, 'Gene_strand'], ';'))
  gene_type_list = unlist(stringr::str_split(sense_Chimera[entry, 'Gene_type'], ';'))
  
  counter = 1
  keep_these_genes = vector()
  keep_these_strands = vector()
  keep_these_types = vector()
  for (gene_strand in gene_strand_list){
   
    if (gene_strand == repeat_strand){
      
      keep_these_genes[counter] = gene_list[counter]
      keep_these_strands[counter] = gene_strand_list[counter]
      keep_these_types[counter] = gene_type_list[counter]
      
    }
    
    else {
      
      keep_these_genes[counter] = ''
      keep_these_strands[counter] = ''
      keep_these_types[counter] = ''
      
    }
    
    counter = counter + 1
     
  }
 
  keep_these_genes = keep_these_genes[keep_these_genes != '']
  sense_Chimera[entry, 'Geneid'] = stringr::str_c(keep_these_genes, collapse = ';')
  
  keep_these_strands = keep_these_strands[keep_these_strands != '']
  sense_Chimera[entry, 'Gene_strand'] = stringr::str_c(keep_these_strands, collapse = ';')
  
  keep_these_types = keep_these_types[keep_these_types != '']
  sense_Chimera[entry, 'Gene_type'] = stringr::str_c(keep_these_types, collapse = ';')
   
}

sense_Chimera = subset(sense_Chimera, Gene_strand != '')
Chimera = sense_Chimera  

################################################################################
# Filtering by expression
################################################################################

## Filter for total occurrence 
Filtered_Chimera = subset(Chimera, total_occurence >= 2) %>%
  mutate(HI_specific = case_when(LO_occurence <= 1 & HI_occurence >= 2 ~ T,
                                 T ~ F)) %>%
  mutate(RepeatID = forcats::fct_reorder(RepeatID, Contribution, .fun='mean'))

## Filter for induction detected by TElocal + DEseq2

#results_df_local_ERE = readRDS(file = glue::glue('{working_directory}/R_variables/results_df_local_ERE'))
#up_TEs = subset(results_df_local_TE, significant == T & log2FoldChange > 0)$locus
#Filtered_Chimera = subset(Chimera, TEid %in% up_TEs) %>%
#  mutate(RepeatID = forcats::fct_reorder(RepeatID, Contribution, .fun='mean'))

## Filter for induction detected by TElocal + edgeR

#edgeR_results_local_ERE = readRDS(file = glue::glue('{working_directory}/R_variables/edgeR_results_local_ERE'))
#up_TEs = subset(edgeR_results_local_ERE, significant == T & logFC > 0)$locus
#Filtered_Chimera = subset(Chimera, TEid %in% up_TEs) %>%
#  mutate(RepeatID = forcats::fct_reorder(RepeatID, Contribution, .fun='mean'))

contribution_plot = ggplot(data = Filtered_Chimera, aes(x = RepeatID, y = Contribution)) +
  geom_bar(stat = 'summary') +
  geom_hline(yintercept = 1, linetype = 'dashed')

contribution_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                     plot.subtitle = element_text(size = 14),
                                     axis.text.x = element_text(size = 10, angle = 90),
                                     axis.text.y = element_text(size = 14),
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_text(size = 15, margin = margin(r = 7.5)),
                                     axis.line = element_line(size = 0.8),
                                     panel.border = element_blank(),
                                     legend.text = element_text(size = 15),
                                     legend.title = element_text(size = 18),
                                     legend.position = c(0.25, 0.93),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/LIONS_contribution_1.png", 
       width = 5, height = 6, units = "in")

################################################################################
#  Contribution analysis: gene-centric
################################################################################

mean_contribution = group_by(Filtered_Chimera, RepeatID) %>% 
  summarize(mean_contribution = mean(Contribution)) %>%
  subset(is.na(RepeatID) == F)
mean_contribution = as.data.frame(mean_contribution)

counter = 1
gene_list = vector()
contribution = vector()
repeat_class = vector()
repeat_name = vector()
LO_occurence = vector()
HI_occurence = vector()
HI_specific = vector()

for (i in 1:nrow(mean_contribution)){
  
  repeat_id = mean_contribution[i, 'RepeatID']
  
  genes = names(sort(table(subset(Filtered_Chimera, RepeatID == repeat_id)$Geneid), decreasing = TRUE))
  
  if (length(genes) > 1){
    
    longest = which.max(stringr::str_length(genes))
    
    genes = genes[longest]
    
  }
  
  specific = names(sort(table(subset(Filtered_Chimera, RepeatID == repeat_id)$HI_specific), decreasing = TRUE))[1]

  genes = unlist(stringr::str_split(genes, pattern = ';'))
  
  for (a in 1:length(genes)){
    gene = stringr::str_trim(genes[a], side = 'both') 
    gene_list[counter] = gene
    HI_specific[counter] = specific
    contribution[counter] = mean_contribution[i, 'mean_contribution']
    repeat_class[counter] = names(sort(table(subset(Filtered_Chimera, RepeatID == repeat_id)$class), decreasing = TRUE)[1])
    LO_occurence[counter] = names(sort(table(subset(Filtered_Chimera, RepeatID == repeat_id)$LO_occurence), decreasing = TRUE)[1])
    HI_occurence[counter] = names(sort(table(subset(Filtered_Chimera, RepeatID == repeat_id)$HI_occurence), decreasing = TRUE)[1])
    repeat_name[counter] = names(sort(table(subset(Filtered_Chimera, RepeatID == repeat_id)$repeatName), decreasing = TRUE)[1])
    
    counter = counter + 1
    
  }
  
}

LIONS_genes = data.frame(ensembl_gene_id = gene_list, contribution = contribution,
                         repeat_name = repeat_name, repeat_class = repeat_class, 
                         LO_occurence = LO_occurence, HI_occurence = HI_occurence,
                         HI_specific = HI_specific) %>%
  mutate(ensembl_gene_id = forcats::fct_reorder(ensembl_gene_id, contribution, .fun='mean'))

## Annotate

ensembl = useEnsembl(biomart = "genes", dataset = 'hsapiens_gene_ensembl')

conversion = getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'external_gene_name'), 
                   mart = ensembl, filter = 'ensembl_gene_id', values = LIONS_genes$ensembl_gene_id)

LIONS_genes_annotated = merge(LIONS_genes, conversion, by = 'ensembl_gene_id', all.x = T) 

LIONS_genes_annotated = mutate(LIONS_genes_annotated, gene_biotype = case_when(gene_biotype == "transcribed_unitary_pseudogene" | gene_biotype == "transcribed_unprocessed_pseudogene" ~ 'Pseudogene',
                                               gene_biotype == 'protein_coding' ~ 'Protein coding',
                                               gene_biotype == 'lncRNA' ~ 'lncRNA',
                                               T ~ gene_biotype))

#LIONS_genes_annotated = mutate(LIONS_genes_annotated, TRA_annotation = case_when(ensembl_gene_id %in% TRA_genes$ensembl_gene_id ~ 'TRA',
#                                      ensembl_gene_id %in% housekeeping_genes$ensembl_gene_id ~ 'Housekeeping',
#                                      T ~ T))

write.csv(LIONS_genes_annotated, '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/gene_lists/LIONS_genes.csv')

## Visualize

# All events
input = LIONS_genes_annotated

# mTEC-HI specific
input = filter(LIONS_genes_annotated, HI_specific == T)

r = "#FBB4AE"
g = "#CCEBC5"
p = "#DECBE4"

contribution_plot = ggplot(data = input, aes(x = ensembl_gene_id, y = contribution, fill = gene_biotype)) +
  geom_bar(stat = 'identity', col = 'black') +
  geom_hline(yintercept = 0.1, linetype = 'dashed') +
  geom_text(aes(label = external_gene_name), angle = 90, nudge_y = 0.3, size = 4) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .1))) +
  ylab('Repeat expression relative to gene') +
  labs(fill= "") +
  scale_fill_manual(labels = c('lncRNA', 'Protein coding', 'Pseudogene'), values = c(r, g, p))

contribution_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                       plot.subtitle = element_text(size = 14),
                                       axis.text.x = element_blank(),
                                       axis.text.y = element_text(size = 14),
                                       axis.title.x = element_blank(),
                                       axis.title.y = element_text(size = 15, margin = margin(r = 7.5)),
                                       axis.line = element_line(size = 0.8),
                                       panel.border = element_blank(),
                                       legend.text = element_text(size = 12),
                                       legend.title = element_text(size = 13),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       legend.position = 'top',
                                       axis.ticks.x=element_blank())

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/LIONS_contribution_genes.png", 
       width = 6, height = 7, units = "in")


## Biotype enrichment

all_genes = getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), mart = ensembl)

all_genes = mutate(all_genes, LIONS = case_when(ensembl_gene_id %in% input$ensembl_gene_id ~ T,
                                                T ~ F))

output = generate_contingency(input = all_genes, 
                     condition_A = list('LIONS' = "LIONS == T",
                                        'Not LIONS' = "LIONS == F"),
                     condition_B = list('lncRNA' = "gene_biotype == 'lncRNA'",
                                        'protein' = "gene_biotype == 'protein_coding'",
                                        'unprocessed_pseudogene' = "gene_biotype == 'transcribed_unprocessed_pseudogene'"),
                     output_type = 'contingency')

chisq.test(output)
chisq.posthoc.test::chisq.posthoc.test(output)

## Gene-set enrichment

library(enrichR)

genes = LIONS_genes$hgnc_symbol
genes = genes[genes != '']

dbs = c("GO_Molecular_Function_2018", 'GO_Biological_Process_2018')

enrichment_results = enrichr(genes, dbs)

plotEnrich(enrichment_results[[1]])

#################################################################
#  Gene enrichment
#################################################################

counter = 1
gene_name = vector()
classification = vector()
gene_group = vector()
for (entry in 1:nrow(LIONS)){
  
  gene_list = unlist(stringr::str_split(LIONS[entry, 'RefID'], pattern = ';'))
  
  for (gene in gene_list){
    
    gene_name[counter] = stringr::str_trim(gene, side = 'both')
    classification[counter] = LIONS[entry, 'classification']
    
    if(gene_name[counter] %in% TRA_genes$external_gene_name){gene_group[counter] = 'TRA'}
    else if(gene_name[counter] %in% housekeeping_genes$external_gene_name){gene_group[counter] = 'housekeeping'}
    else{gene_group[counter] = 'Other'}
      
    counter = counter + 1
      
  }
    
}
  
LIONS_genes = data.frame(gene_name = gene_name, classification = classification, gene_group = gene_group)

gene_group = c('TRA', 'housekeeping')
p_value = vector()
odds_ratio = vector()
lower_interval = vector()
upper_interval = vector()

for (i in 1:length(gene_group)){
  
  query_group = gene_group[i]
  
  column_1 = c(nrow(subset(LIONS_genes, gene_group == query_group & classification == 'HI')),
               nrow(subset(LIONS_genes, gene_group == query_group & classification == 'LO')))
  
  column_2 = c(nrow(subset(LIONS_genes, gene_group != query_group & classification == 'HI')),
               nrow(subset(LIONS_genes, gene_group != query_group & classification == 'LO')))
  
  contingency = data.frame(query_group = column_1, 'Not' = column_2, row.names = c('HI', 'LO'))
  
  print(query_group)
  print(contingency)
  
  p_value[i] = fisher.test(x = contingency)[[1]]
  odds_ratio[i] = fisher.test(x = contingency)[[3]]
  lower_interval[i] = fisher.test(x = contingency)$conf.int[1]
  upper_interval[i] = fisher.test(x = contingency)$conf.int[2]
  
}

p_value = p.adjust(p_value, method = 'bonferroni')

output = data.frame(gene_group = gene_group, 
                    p_value = p_value, 
                    odds_ratio = odds_ratio, 
                    lower_interval = lower_interval, 
                    upper_interval = upper_interval) %>%
  mutate(significant = case_when(p_value < 0.001 ~ '***', 
                                 p_value < 0.01 ~ '**',
                                 p_value < 0.05 ~ '*',
                                 p_value >= 0.05 ~ '')) %>%
  mutate(class = forcats::fct_reorder(gene_group, odds_ratio))

odds_ratio_plot = ggplot(data = output, aes(x = gene_group, y = odds_ratio)) + 
  geom_point() +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_errorbar(aes(ymin = lower_interval, ymax = upper_interval), width = 0) +
  xlab('') +
  ylab('Odds ratio') +
  geom_text(aes(label = significant), nudge_y = 6, size = 6)
  
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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/LIONS_odds_ratio.png", 
       width = 4, height = 5.25, units = "in")