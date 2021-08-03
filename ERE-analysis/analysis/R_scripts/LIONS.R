library(ggvenn)
library(tidyverse)
library(GenomicRanges)
library(biomaRt)

#################################################################
#  Data import
#################################################################

LIONS_HI = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/External_packages/LIONS/mTEC-analysis/mTEC-analysis.rslions', 
         sep = '\t', header = T)
LIONS_HI['LO_occurence'] = LIONS_HI['Normal_Occ']
LIONS_HI['HI_occurence'] = LIONS_HI['Cancer_Occ']
LIONS_HI = dplyr::select(LIONS_HI, -c('Normal_Occ', 'Cancer_Occ')) %>%
  mutate(classification = 'HI')

LIONS_LO = read.csv(file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/External_packages/LIONS/mTEC-analysis/mTEC-analysis.inv.rslions', 
                    sep = '\t', header = T)
LIONS_LO['LO_occurence'] = LIONS_LO['Cancer_Occ']
LIONS_LO['HI_occurence'] = LIONS_LO['Normal_Occ']
LIONS_LO = dplyr::select(LIONS_LO, -c('Normal_Occ', 'Cancer_Occ')) %>%
  mutate(classification = 'LO')

LIONS = bind_rows(LIONS_HI, LIONS_LO) %>%
  separate(coordinates, into = c('chromosome', 'range'), sep = ':') %>%
  separate(range, into = c('start', 'end'), sep = '-') %>%
  separate(repeatName, remove = F, into = c('sub-family', 'class', 'family'), sep = ':')

#################################################################
#  Contribution analysis
#################################################################

home_directory = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/External_packages/LIONS/mTEC-analysis'

STDIN = c(glue::glue('{home_directory}/mTEC-analysis.lions'),
          glue::glue('{home_directory}/input.list'),
          glue::glue('{home_directory}/rmStats.Rdata'),
          1,
          0)

pLIONS = STDIN[1] # Input LIONS project file
GROUP_LIST = STDIN[2]
RMDATA_PATH= STDIN[3]
MultiCan = as.numeric(STDIN[4])
NonNormal = as.numeric(STDIN[5])

Chimera = read.csv(pLIONS, header=T, sep='\t')
ID = paste(Chimera[,'RepeatID'],Chimera[,'LIBRARY'],sep='~')
IDdup = !duplicated(ID)
Chimera = Chimera[IDdup,]

Chimera_annotated = mutate(Chimera, EStrand = case_when(EStrand == 1 ~ '+',
                                                        EStrand == -1 ~ '-',
                                                        EStrand == 0 ~ '*')) %>%
  mutate(Chromosome = glue::glue('chr{Chromosome}'))

Chimera_annotated = makeGRangesFromDataFrame(Chimera_annotated,
                                             start.field = 'EStart',
                                             end.field = 'EEnd',
                                             seqnames.field = 'Chromosome',
                                             strand.field = 'EStrand',
                                             keep.extra.columns = T)

gene_overlaps = findOverlaps(query = Chimera_annotated,
                             subject = GRanges_gene, 
                             ignore.strand = T)

gene_overlaps = as.data.frame(gene_overlaps)

Geneid = vector()
for (entry in 1:length(Chimera_annotated)){
  
  if (entry %in% gene_overlaps$queryHits){
    
    gene_hits = subset(gene_overlaps, queryHits == entry)[, 'subjectHits']
    gene_hits = GRanges_gene[gene_hits, ]$Geneid
    
    Geneid[entry] = paste(gene_hits, collapse = ';')
    
  }
  
  else{
    
    Geneid[entry] = NA
    
  }
  
}

Chimera['Geneid'] = Geneid

Chimera['Group'] = sapply(Chimera[,'LIBRARY'], FUN = GroupAssign)
Chimera['LO_Occ'] = apply(X = Chimera, MARGIN = 1, FUN = MatchEntry, GroupN = 1)
Chimera['HI_Occ'] = apply(X = Chimera, MARGIN = 1, FUN = MatchEntry, GroupN = 2)
Chimera = mutate(Chimera, total_Occ = LO_Occ + HI_Occ)

## Annotate

TRA_annotation = vector()
gene_type_annotation = vector()
for (entry in 1:nrow(Chimera)){
  
  gene_list = unlist(stringr::str_split(Chimera[entry, 'Geneid'], pattern = ';'))

  for (gene in 1:length(gene_list)){
    
    gene_entry = stringr::str_trim(gene_list[gene], side = 'both')
    gene_list[gene] = gene_entry
    
  }
  
  if (T %in% (gene_list %in% TRA_genes$ensembl_gene_id)){
    
    TRA_annotation[entry] = T
    
  }
  
  else{
    
    TRA_annotation[entry] = F
    
  }
  
  
}

Chimera['TRA_annotation'] = TRA_annotation

## Strict filtering
Filtered_Chimera = subset(Chimera, Group == 2) %>%
  subset(LO_Occ <= 1 & HI_Occ >= 2) %>%
  mutate(RepeatID = forcats::fct_reorder(RepeatID, Contribution, .fun='mean'))

## Not-so-strict 
Filtered_Chimera = subset(Chimera, Group == 2) %>%
  subset(LO_Occ <= 0 & HI_Occ >= 1) %>%
  mutate(RepeatID = forcats::fct_reorder(RepeatID, Contribution, .fun='mean'))

contribution_plot = ggplot(data = Filtered_Chimera, aes(x = RepeatID, y = Contribution, fill = TRA_annotation)) +
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

#################################################################
#  Gene enrichment
#################################################################

mean_contribution = group_by(Filtered_Chimera, RepeatID) %>% 
  summarize(mean_contribution = mean(Contribution))
mean_contribution = as.data.frame(mean_contribution)

counter = 1
gene_list = vector()
contribution = vector()
for (i in 1:nrow(mean_contribution)){
  
  repeat_id = mean_contribution[i, 'RepeatID']
  
  genes = subset(Filtered_Chimera, RepeatID == repeat_id)[1, 'Geneid']
  
  genes = unlist(stringr::str_split(genes, pattern = ';'))
  print(genes)
  
  for (gene in genes){
    
    gene_list[counter] = gene
    contribution[counter] = mean_contribution[i, 'mean_contribution']
    counter = counter + 1
    
  }
  
}

LIONS_genes = data.frame(ensembl_gene_id = gene_list, contribution = contribution) %>%
  subset(!is.na(ensembl_gene_id)) %>%
  mutate(ensembl_gene_id = forcats::fct_reorder(ensembl_gene_id, contribution, .fun='mean'))

conversion = getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), filters = 'ensembl_gene_id', values = LIONS_genes$ensembl_gene_id, mart= ensembl)

LIONS_genes = merge(LIONS_genes, conversion , by = 'ensembl_gene_id') 

## Annotate

LIONS_genes = mutate(LIONS_genes, TRA_annotation = case_when(ensembl_gene_id %in% TRA_genes$ensembl_gene_id ~ T,
                                           !(ensembl_gene_id %in% TRA_genes$ensembl_gene_id) ~ F))

gene_type = vector()
for (i in 1:nrow(LIONS_genes)){
  
  gene = LIONS_genes[i, 'ensembl_gene_id']
  gene_type[i] = subset(GRanges_gene, Geneid == gene)$Class
  
}

LIONS_genes['gene_type'] = gene_type

contribution_plot = ggplot(data = LIONS_genes, aes(x = ensembl_gene_id, y = contribution, fill = TRA_annotation)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = .1, linetype = 'dashed') +
  geom_text(aes(label = hgnc_symbol, angle = 270), nudge_y = 0.4)

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

ggsave("/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/LIONS_contribution_genes.png", 
       width = 9, height = 7, units = "in")

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

#################################################################
#  ?
#################################################################

LIONS_overlaps = findOverlaps(query = GRanges_LIONS, subject = GRanges_ERE)
LIONS_overlaps = as.data.frame(LIONS_overlaps)

for (i in 1:nrow(LIONS)){
  
  ERE = LIONS[i, 'sub-family']
  print(ERE)
  
}

## Venn diagram

venn_diagram_list = list('mTEC-LO' = subset(LIONS, Normal_Occ >= 2)$transcriptID,
                         'mTEC-HI' = subset(LIONS, Cancer_Occ >= 2)$transcriptID)

ggvenn(venn_diagram_list)

## Genes

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
acceptable_genes = getBM(attributes=c("hgnc_symbol"), mart= ensembl)$hgnc_symbol

gene_name = vector()
classification = vector()
counter = 1
for (entry in 1:nrow(LIONS)){
  
  gene_list = unlist(stringr::str_split(LIONS[entry, 'RefID'], pattern = ';'))
  
  for (gene in gene_list){
    
    gene_entry = stringr::str_trim(gene, side = 'both')
    
    if (gene_entry %in% acceptable_genes){
      
      gene_name[counter] = gene_entry
      
      if (LIONS[entry, 'Normal_Occ'] == 0 & LIONS[entry, 'Cancer_Occ'] >= 1){
        
        classification[counter] = 'HI'
        
      }
      
      else if (LIONS[entry, 'Normal_Occ'] >= 1 & LIONS[entry, 'Cancer_Occ'] == 0){
        
        classification[counter] = 'LO'
        
      }
      
      else{
        
        classification[counter] = 'NONE'
        
      }
      
      counter = counter + 1
      
    }

  }
  
}

LIONS_genes = data.frame(hgnc_symbol = gene_name, classification = classification) %>%
  distinct() 

conversion = getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), filters = 'hgnc_symbol', values = LIONS_genes$hgnc_symbol, mart= ensembl)

LIONS_genes = merge(LIONS_genes, conversion , by = 'hgnc_symbol') %>%
  mutate(gene_group = case_when(ensembl_gene_id %in% AIRE_genes ~ 'AIRE',
                               ensembl_gene_id %in% housekeeping_genes ~ 'Housekeeping',
                               ensembl_gene_id %in% FEZF2_genes ~ 'FEZF2',
                               !(ensembl_gene_id %in% AIRE_genes) & !(ensembl_gene_id %in% housekeeping_genes) & !(ensembl_gene_id %in% FEZF2_genes) ~ 'Other'))

