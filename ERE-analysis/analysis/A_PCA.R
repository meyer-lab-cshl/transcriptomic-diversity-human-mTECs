library(DESeq2)
library(ggplot2)

extract_subset = function(mode, input){
  
  ##############################################################################
  # Takes a dataframe or DeSeq2 object ('input') and extracts a subset of elements
  # depending on 'mode'. 'gene' extracts genes, 'TE' extracts all annotated repeats,
  # and 'ERE' extracts endogenous retroelements (i.e. LTRs, LINEs, SINEs)
  ##############################################################################
  
  if (mode == 'gene'){
    
    output = input[grep("^ENSG",rownames(input)),]
    
  }
  
  if (mode == 'TE'){
    
    output = input[grepl("^(?!ENSG).*$",rownames(input), perl = TRUE),]
    
  }
  
  if (mode == 'ERE'){
    
    output = input[grepl("^(?!ENSG).*$",rownames(input), perl = TRUE),]
    output = output[!grepl("Satellite",rownames(output), perl = TRUE),]
    output = output[!grepl("DNA",rownames(output), perl = TRUE),]
    output = output[!grepl("DNA",rownames(output), perl = TRUE),]
    
  }
  
  return(output)
  
}

differential_expression = function(raw_count_table, min_reads = 2, design = ~tissue){
  
  ##############################################################################
  # Generates a DESeq2 dds object.  By default, design 
  # detects differences between tissues. For differential expression analysis 
  # between mTEC-HI and -LO samples, you'll want to change the design to
  # '~patient + tissue', but this can't be done when including GTEx data. After
  # producing the dds object, filters lowly expressed genes/TEs by requiring
  # a minimum number of reads (default 2) in a row.
  ##############################################################################
  
  ## Define sampleInfo
  
  ID = colnames(raw_count_table)
  sampleInfo = data.frame(ID,row.names=colnames(raw_count_table))
  sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'tissue', 'batch'), sep = '_'))
  sampleInfo$patient = factor(sampleInfo$patient)
  
  ## Construct DESeq dataset object
  
  dds <- DESeqDataSetFromMatrix(countData = raw_count_table, colData = sampleInfo, design = design)
  
  ## Run differential expression analysis
  
  dds = DESeq(dds)
  
  ## Pre-filter to remove rows with less than 'min_reads' total (default 2)
  
  dds = dds[ rowSums(counts(dds)) >= min_reads, ]
  
  return(dds)
  
}

############################################################################################
##                                        PCA 
############################################################################################

## 'count_table_directory' should contain the text file 'TE_transcripts_counts' containing
## the raw counts from each tissue of interest. Each entry should be labelled in the format 
## '{unique ID}_{tissue}_{batch}'. e.g. 'pt214_mTEC-hi_our-data'

## Data import:

count_table_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/"
data = read.table(glue('{count_table_directory}TE_transcripts_counts'),header=T,row.names=1)

## Run DESeq2 to obtain normalized counts

TE_data = extract_subset(mode = 'TE', input = data)

dds_transcripts_TE = differential_expression(TE_data, design=~tissue)

vs_dds_transcripts_TE = vst(dds_transcripts_TE, blind=FALSE)

assay(vs_dds_transcripts_TE) = limma::removeBatchEffect(assay(vs_dds_transcripts_TE), vs_dds_transcripts_TE$batch)

## Plot

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