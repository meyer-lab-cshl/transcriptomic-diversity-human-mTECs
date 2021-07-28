library(DESeq2)
library(tidyverse)
library(ggplot2)
library(glue)

functions_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_functions/"
functions = c('extract_subset', 'differential_expression', 'process_DESeq2_results', 'build_count_table', 'make_GRanges')

for (i in functions){
  
  load(glue('{functions_directory}{i}'))
  
}

#################################################################
# Salmon
#################################################################

files = list.files(path="~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/SalmonTE", 
                   pattern="*raw.csv", 
                   full.names=TRUE, 
                   recursive=FALSE)

counter = 0
for (file in files){
  
  input = read.csv(file = file)
  if (counter == 0){mTEC_data = input}
  else{mTEC_data = merge(mTEC_data, input, by = 'TE')}
  counter = counter + 1
  
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

#################################################################
# TEtranscripts
#################################################################

## 'count_table_directory' should contain the text file 'TE_transcripts_hi_vs_lo.cntTable' containing
## the raw counts for mTEC-HI and -LO conditions. Each entry should be labelled in the format 
## '{unique ID}_{tissue}_{batch}'. e.g. 'pt214_mTEC-hi_our-data'

count_table_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/"
data = read.table(glue('{count_table_directory}mTEC.cntTable'),header=T,row.names=1)

mTEC_data = read.table(glue('{count_table_directory}TE_transcripts_counts'),header=T,row.names=1) %>%
  select(c('pt214_mTEC.lo_new_Aligned.out.bam', 'pt226_mTEC.lo_new_Aligned.out.bam', 'pt221_mTEC.lo_new_Aligned.out.bam',
           'pt214_mTEC.hi_new_Aligned.out.bam', 'pt226_mTEC.hi_new_Aligned.out.bam', 'pt221_mTEC.hi_new_Aligned.out.bam'))

## Run DESeq2

ERE_data = extract_subset(mTEC_data, mode = 'ERE')

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
        file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/results_df_transcripts_ERE')

#keep_list = c('LTR', 'LINE', 'SINE')

#results_df_transcripts_TE = process_DESeq2_results(results = results_transcripts_TE, mode = 'TE_transcripts') %>%
#  mutate(ID = sub("\\?", "", ID)) %>%
#  mutate(class = sub("\\?", "", class)) %>%
#  mutate(class = case_when(class %in% keep_list ~ class,
#                           !(class %in% keep_list) ~ 'Other'))


#################################################################
# TElocal
#################################################################

count_table_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/"
data = read.table(glue::glue('{count_table_directory}TE_local_hi_vs_lo.cntTable'),header=T,row.names=1)

## DESeq2

dds_local = differential_expression(data, design=~patient+tissue)

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
        file = '~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_variables/results_df_local_ERE')

## TPM-like counts

hi_raw_counts = select(data, c('pt226_hi_fastp_1.fastq_Aligned.out.bam', 'pt214_hi_fastp_1.fastq_Aligned.out.bam', 'pt221_hi_fastp_1.fastq_Aligned.out.bam'))

hi_ERE_raw_counts = extract_subset(input = hi_raw_counts, mode = 'ERE')
hi_ERE_raw_counts = cbind(ID = rownames(hi_ERE_raw_counts), hi_ERE_raw_counts) %>%
  tidyr::separate(col = 'ID', into = c('locus', 'gene', 'family', 'class'), sep = ':')
hi_ERE_raw_counts$ID = cbind(ID = rownames(hi_ERE_raw_counts), hi_ERE_raw_counts) 

hi_ERE_raw_counts_annotated = make_GRanges(results_df = hi_ERE_raw_counts, mode = 'TE')
hi_ERE_raw_counts_annotated$width = width(hi_ERE_raw_counts_annotated)
hi_ERE_raw_counts_annotated = as.data.frame(hi_ERE_raw_counts_annotated)

estimate_TPM = function(input, sample_names){
  
  for (name in sample_names){
    
    print(name)
    
    input = mutate(input, name = name / width)
    
    colsum = sum(input$name)
    
    input = mutate(input, name = name * 1e6 / colsum)
    
  }
  
  return(input)

}
estimate_TPM(hi_ERE_raw_counts_annotated, sample_names = c('pt226_hi_fastp_1.fastq_Aligned.out.bam',
                                                           'pt214_hi_fastp_1.fastq_Aligned.out.bam',
                                                           'pt221_hi_fastp_1.fastq_Aligned.out.bam'))

input = hi_ERE_raw_counts_annotated

input = mutate(input, pt226_hi_fastp_1.fastq_Aligned.out.bam = pt226_hi_fastp_1.fastq_Aligned.out.bam / width)
colsum = sum(input$pt226_hi_fastp_1.fastq_Aligned.out.bam)
input = mutate(input, pt226_hi_fastp_1.fastq_Aligned.out.bam = pt226_hi_fastp_1.fastq_Aligned.out.bam * 1e6 / colsum)

input = mutate(input, pt214_hi_fastp_1.fastq_Aligned.out.bam = pt214_hi_fastp_1.fastq_Aligned.out.bam / width)
colsum = sum(input$pt214_hi_fastp_1.fastq_Aligned.out.bam)
input = mutate(input, pt214_hi_fastp_1.fastq_Aligned.out.bam = pt214_hi_fastp_1.fastq_Aligned.out.bam * 1e6 / colsum)

input = mutate(input, pt221_hi_fastp_1.fastq_Aligned.out.bam = pt221_hi_fastp_1.fastq_Aligned.out.bam / width)
colsum = sum(input$pt221_hi_fastp_1.fastq_Aligned.out.bam)
input = mutate(input, pt221_hi_fastp_1.fastq_Aligned.out.bam = pt221_hi_fastp_1.fastq_Aligned.out.bam * 1e6 / colsum)

hi_ERE_TPM = input

matrix = select(hi_ERE_TPM, c('pt226_hi_fastp_1.fastq_Aligned.out.bam',
                              'pt214_hi_fastp_1.fastq_Aligned.out.bam',
                              'pt221_hi_fastp_1.fastq_Aligned.out.bam'))

hi_ERE_TPM$mean_TPM = rowMeans(matrix)

hi_ERE_TPM_filtered = subset(hi_ERE_TPM, mean_TPM >= 10)

hi_ERE_TPM_filtered = mutate(hi_ERE_TPM_filtered, ereMAP = case_when(ID.ID %in% ereMAP_loci ~ T,
                                                                     !(ID.ID %in% ereMAP_loci) ~ F)) %>%
  mutate(locus = forcats::fct_reorder(locus, mean_TPM))

plot =ggplot(data = subset(hi_ERE_TPM_filtered, ereMAP == T), aes(x = locus, y = log10(mean_TPM))) +
  geom_point(aes(fill = ereMAP), shape = 21, stroke = 0, size = 2)

plot + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                                  plot.subtitle = element_text(size = 14),
                                  axis.text.x = element_blank(),
                                  axis.text.y = element_text(size = 14),
                                  axis.title = element_text(size = 14),
                                  axis.line = element_line(size = 0.8),
                                  panel.border = element_blank(),
                                  legend.text = element_text(size = 15),
                                  legend.title = element_text(size = 18),
                                  legend.position = c(0.2, 0.93),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())

subset(hi_ERE_TPM, ID.ID %in% ereMAP_loci)

## Test only ereMAPs

ereMAP_data = data[row.names(data) %in% ereMAP_loci, ]

dds_local_ereMAP = differential_expression(ereMAP_data, design=~patient+tissue)

results_local_ereMAP = results(dds_local_ereMAP, 
                        independentFiltering = F,
                        contrast = c('tissue', 'hi', 'lo'))

results_df_local_ereMAP = process_DESeq2_results(results = results_local_ereMAP, mode = 'TE_local')

