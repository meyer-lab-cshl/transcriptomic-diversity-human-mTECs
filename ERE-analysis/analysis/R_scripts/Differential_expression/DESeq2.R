################################################################################
# DESeq2 for differential expression analysis of TE_count/TE_local/SalmonTE 
# output
################################################################################

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(glue)

working_directory = '/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis'
functions_directory = glue("{working_directory}/R_functions/")

functions = c('extract_subset', 'differential_expression', 'process_DESeq2_results', 'build_count_table', 'make_GRanges')

for (i in functions){
  
  load(glue('{functions_directory}{i}'))
  
}

################################################################################
# SalmonTE
################################################################################

files = list.files(path=glue("{working_directory}/count_tables/SalmonTE"), 
                   pattern="*mTEC*", 
                   full.names=TRUE, 
                   recursive=FALSE)

counter = 0
for (file in files){
  
  input = read.csv(file = file)
  if (counter == 0){mTEC_data = input}
  else{mTEC_data = merge(mTEC_data, input, by = 'TE')}
  counter = counter + 1
  remove(input)
  
}

row.names(mTEC_data) = mTEC_data$TE
mTEC_data = select(mTEC_data, -TE)

mTEC_data_rounded = round(mTEC_data)

Salmon_dds_transcripts = differential_expression(mTEC_data_rounded, 
                                                 design=~patient + tissue)

Salmon_vs_dds_transcripts = vst(Salmon_dds_transcripts, blind=FALSE)

Salmon_results_transcripts = results(Salmon_dds_transcripts, 
                              contrast = c('tissue', 'mTEC.hi', 'mTEC.lo'), 
                              independentFiltering = F)

Salmon_results_transcripts_df = process_DESeq2_results(results = Salmon_results_transcripts, mode = 'Salmon')
Salmon_results_transcripts_df$TE = row.names(Salmon_results_transcripts_df)

################################################################################
# TEtranscripts
################################################################################

## 'count_table_directory' should contain a text file with the raw counts for 
## mTEC-HI and -LO conditions. Each entry should be labelled in the format 
## '{unique ID}_{tissue}_{batch}'. e.g. 'pt214_mTEC-hi_our-data'

count_table_directory = glue("{working_directory}/count_tables/TE_transcripts")
data = read.table(glue('{count_table_directory}/mTEC.cntTable'),header=T,row.names=1)

ERE_data = extract_subset(mTEC_data, mode = 'ERE')

## Run DESeq2

dds_transcripts = differential_expression(ERE_data, design=~patient + tissue)
vs_dds_transcripts = vst(dds_transcripts, blind=FALSE)

results_transcripts = results(dds_transcripts, 
                              contrast = c('tissue', 'mTEC.hi', 'mTEC.lo'), 
                              independentFiltering = F)

results_df_transcripts_ERE = process_DESeq2_results(results = results_transcripts_ERE, 
                                                    mode = 'TE_transcripts') %>%
  mutate(ID = sub("\\?", "", ID)) %>%
  mutate(class = sub("\\?", "", class))

saveRDS(object = results_df_transcripts_ERE, 
        file = glue('{working_directory}/R_variables/results_df_transcripts_ERE'))

################################################################################
# TElocal
################################################################################

count_table_directory = glue("{working_directory}/count_tables/TE_local")
data = read.table(glue::glue('{count_table_directory}/TE_local_hi_vs_lo.cntTable'),header=T,row.names=1)

## DESeq2

dds_local = differential_expression(data, design=~patient+tissue, min_reads = 2)

vs_dds_local = vst(dds_local, blind = F)

results_local = results(dds_local, 
                        independentFiltering = F,
                        contrast = c('tissue', 'hi', 'lo'))

results_local_gene = extract_subset(mode = 'gene', input = results_local)
results_local_TE = extract_subset(mode = 'TE', input = results_local)
results_local_ERE = extract_subset(mode = 'ERE', input = results_local)

results_df_local_gene = process_DESeq2_results(results = results_local_gene, mode = 'Gene') 
results_df_local_gene$Geneid = gsub('\\..+$', '', results_df_local_gene$Geneid)

keep_list = c('LTR', 'LINE', 'SINE', 'Satellite', 'DNA')

results_df_local_TE = process_DESeq2_results(results = results_local_TE, mode = 'TE_local') %>%
  mutate(ID = sub("\\?", "", ID)) %>%
  mutate(class = sub("\\?", "", class)) %>%
  mutate(class = case_when(class %in% keep_list ~ class,
                           !(class %in% keep_list) ~ 'Other'))

results_df_local_ERE = process_DESeq2_results(results = results_local_ERE, mode = 'TE_local') %>%
  mutate(ID = sub("\\?", "", ID)) %>%
  mutate(class = sub("\\?", "", class))

saveRDS(object = results_df_local_ERE, 
        file = glue('{working_directory}/R_variables/results_df_local_ERE'))

## Test only ereMAPs

ereMAP_data = data[row.names(data) %in% ereMAP_loci, ]

dds_local_ereMAP = differential_expression(ereMAP_data, design=~patient+tissue)

results_local_ereMAP = results(dds_local_ereMAP, 
                        independentFiltering = F,
                        contrast = c('tissue', 'hi', 'lo'))

results_df_local_ereMAP = process_DESeq2_results(results = results_local_ereMAP, mode = 'TE_local')

saveRDS(object = results_df_local_ereMAP, 
        file = glue('{working_directory}/R_variables/results_df_local_ereMAP'))

