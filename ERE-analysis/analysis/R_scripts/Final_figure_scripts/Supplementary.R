library(tidyverse)
library(glue)
library(pheatmap)

working_directory = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/'
functions_directory = glue("{working_directory}/R_functions")
functions = c('generate_contingency', 'generate_heatmap_matrix')

for (i in functions){
  
  load(glue('{functions_directory}/{i}'))
  
}

################################################################################
# Enrichment of LTRs (supplementary A)
################################################################################

## Import required variables
results_df_transcripts_ERE = readRDS(file = glue('{working_directory}/R_variables/results_df_transcripts_ERE')) 

## Chi-square test

output = generate_contingency(input = results_df_transcripts_ERE,
                     condition_A = list('LTR' = "class == 'LTR'", 
                                     'SINE' = "class == 'SINE'", 
                                     'LINE' = "class == 'LINE'",
                                     'Retroposon' = "class == 'Retroposon'"),
                     condition_B = list('up' = 'significant == T & log2FoldChange > 0',
                                        'unchanged' = 'significant == F',
                                        'down' = 'significant == T & log2FoldChange < 0 '),
                     output_type = 'contingency')

chisq.test(output)
chisq.posthoc.test::chisq.posthoc.test(output)

## Bar chart

r = "#e41a1c"
b = '#377eb8'
g = "#4daf4a"
p = "#984ea3"

output = generate_contingency(input = results_df_transcripts_ERE,
                              condition_A = list('LTR' = "class == 'LTR'", 
                                                 'SINE' = "class == 'SINE'", 
                                                 'LINE' = "class == 'LINE'",
                                                 'Retroposon' = "class == 'Retroposon'"),
                              condition_B = list('up' = 'significant == T & log2FoldChange > 0',
                                                 'unchanged' = 'significant == F',
                                                 'down' = 'significant == T & log2FoldChange < 0 '),
                              output_type = 'count')

output$condition_A = factor(output$condition_A, levels = c('LTR', 'SINE', 'LINE', 'Retroposon'))

bar_chart = ggplot(output, aes(x = condition_B, y = percent, fill = condition_A)) + 
  geom_col(colour = 'black', position = 'fill') +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  xlab('EREs') +
  ylab('Fraction') +
  labs(fill= "") +
  scale_x_discrete(labels = c('Down', 'Unchanged', 'Up')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_fill_manual(labels = c('LTR', 'SINE', 'LINE', 'Composite'), values = c(g, p, r, b))

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

################################################################################
# ereMAP expression (supplementary B)
################################################################################

LTR = subset(counts_annotated, class == 'LTR')$locus
LINE = subset(counts_annotated, class == 'LINE')$locus
SINE = subset(counts_annotated, class == 'SINE')$locus

## Expression above RPKM thresholds

input = readRDS(file = glue('{working_directory}/R_variables/ereMAP_mean_RPKM'))

thresholds = c(0.1, 1, 10, 1000)
counter_1 = 0
df = data.frame(row.names = colnames(input))
for (threshold in thresholds){
  
  counter_1 = counter_1 + 1
  RPKM_threshold = vector()
  counter_2 = 0
  for (tissue in colnames(input)){
    
    counter_2 = counter_2 + 1
    RPKM_threshold[counter_2] = nrow(subset(input, eval(parse(text=tissue)) > threshold))
    
  }
  
  df[, counter_1] = RPKM_threshold
  
}
names(df) = thresholds


## Boxplot of RPKM values

input = readRDS(file = glue('{working_directory}/R_variables/ereMAP_mean_RPKM'))
for (name in 1:length(colnames(input))){
  
  if(colnames(input)[name] == 'Small.Intestine.Terminal.Ileum'){
    
    colnames(input)[name] = 'Small intestine'
    
  }
  
  if(colnames(input)[name] == 'Brain.Caudate'){
    
    colnames(input)[name] = 'Basal ganglia'
    
  }
  
  if(colnames(input)[name] == 'Brain.Frontal.Cortex'){
    
    colnames(input)[name] = 'Frontal cortex'
    
  }
  
  
  if(colnames(input)[name] == 'Adrenal.Gland'){
    
    colnames(input)[name] = 'Adrenal gland'
    
  }
  
  if(colnames(input)[name] == 'Brain.Cerebellum'){
    
    colnames(input)[name] = 'Cerebellum'
    
  }
  
  if(colnames(input)[name] == 'ESC'){
    
    colnames(input)[name] = 'ESC'
    
  }
  
  if(colnames(input)[name] == 'mTEC.hi'){
    
    colnames(input)[name] = 'mTEC-hi'
    
  }
  
  if(colnames(input)[name] == 'mTEC.lo'){
    
    colnames(input)[name] = 'mTEC-lo'
    
  }
  
  if(colnames(input)[name] == 'Colon.Transverse'){
    
    colnames(input)[name] = "Colon"
    
  }
  
  if(colnames(input)[name] == 'Brain.Spinal.Cord'){
    
    colnames(input)[name] = "Spinal cord"
    
  }
  
  if(colnames(input)[name] == 'Heart.Left.Ventricle'){
    
    colnames(input)[name] = "Heart"
    
  }
  
  if(colnames(input)[name] == 'Kidney.Cortex'){
    
    colnames(input)[name] = "Kidney"
    
  }
  
  if(colnames(input)[name] == 'Muscle.Skeletal'){
    
    colnames(input)[name] = "Skeletal muscle"
    
  }
  
  if(colnames(input)[name] == 'Breast.Mammary.Tissue'){
    
    colnames(input)[name] = "Breast"
    
  }
  
  if(colnames(input)[name] == 'Esophagus.Mucosa'){
    
    colnames(input)[name] = "Esophagus"
    
  }
  
  if(colnames(input)[name] == 'Adipose.Subcutaneous'){
    
    colnames(input)[name] = "Adipose"
    
  }
  
  if(colnames(input)[name] == 'Brain.Substantia.nigra'){
    
    colnames(input)[name] = "Substantia nigra"
    
  }
  
  if(colnames(input)[name] == 'Skin.Not.Sun.Exposed'){
    
    colnames(input)[name] = "Skin (unexposed)"
    
  }
  
  if(colnames(input)[name] == 'Skin.Sun.Exposed'){
    
    colnames(input)[name] = "Skin (sun exposed)"
    
  }
  
}

input$locus = rownames(input)
input = mutate(input, class = case_when(locus %in% LTR ~ 'LTR',
                                        locus %in% LINE ~ 'LINE',
                                        locus %in% SINE ~ 'SINE'))

input = pivot_longer(input, cols = c(1:28), names_to = 'tissue')

input = mutate(input, color = case_when(tissue == 'mTEC-hi' ~ '#4c72b0ff',
                                 tissue == 'mTEC-lo' ~ '#dd8452ff',
                                 T ~ 'white'))

box_plot = ggplot(data = input, 
                  aes(x = reorder(tissue,value, FUN = median), 
                      y = log2(value+0.01), fill = color)) +
  geom_boxplot() +
  xlab('Tissue') +
  ylab(expression('log'[2]*'(RPKM)')) +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff', '#D9D9D9'))

box_plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                               plot.subtitle = element_text(size = 14),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text.x = element_text(size = 13, 
                                                          angle = 90,
                                                          hjust=0.97, vjust=0.2),
                               axis.text.y = element_text(size = 14),
                               axis.title.y = element_text(size = 14),
                               axis.title.x = element_text(size = 14, margin = margin(t = 6)),
                               axis.line = element_line(size = 0.8),
                               panel.border = element_blank(),
                               legend.position = 'none')

## Heatmap

input = readRDS(file = glue('{working_directory}/R_variables/ereMAP_mean_RPKM'))

my_heatmap = pheatmap(input, 
                      cluster_rows=T,
                      show_rownames=F,
                      show_colnames = T,
                      cluster_cols=T,
                      scale = 'row',
                      angle_col = 45)

