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
#  Gene enrichment
#################################################################

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
acceptable_genes = getBM(attributes=c("hgnc_symbol"), mart= ensembl)$hgnc_symbol

counter = 1
gene_name = vector()
classification = vector()
for (entry in 1:nrow(LIONS)){
  
  gene_list = unlist(stringr::str_split(LIONS[entry, 'RefID'], pattern = ';'))
  
  for (gene in gene_list){
    
    gene_entry = stringr::str_trim(gene, side = 'both')
    
    if (gene_entry %in% acceptable_genes){
      
      gene_name[counter] = gene_entry
      classification[counter] = LIONS[entry, 'classification']
      
      counter = counter + 1
      
    }
    
  }
  
}

LIONS_genes = data.frame(hgnc_symbol = gene_name, classification = classification)

conversion = getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), filters = 'hgnc_symbol', values = LIONS_genes$hgnc_symbol, mart= ensembl)

LIONS_genes = merge(LIONS_genes, conversion, by = 'hgnc_symbol') %>%
  mutate(gene_group = case_when(ensembl_gene_id %in% TRA_genes ~ 'TRA',
                                ensembl_gene_id %in% housekeeping_genes ~ 'Housekeeping',
                                !(ensembl_gene_id %in% TRA_genes) & !(ensembl_gene_id %in% housekeeping_genes) ~ 'Other')) %>%
  mutate(gene_expression_group = case_when(ensembl_gene_id %in% up_genes ~ 'Up',
                                           ensembl_gene_id %in% down_genes ~ 'Down',
                                           !(ensembl_gene_id %in% down_genes) & !(ensembl_gene_id %in% up_genes) ~ 'Unchanged'))

gene_group = c('TRA', 'Housekeeping')
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

