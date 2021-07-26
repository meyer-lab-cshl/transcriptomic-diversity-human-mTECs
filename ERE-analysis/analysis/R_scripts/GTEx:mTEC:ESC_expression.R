library(pheatmap)
library(tidyverse)

functions_directory = "/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/R_functions/"
functions = c('save_pheatmap_png')

for (i in functions){
  
  load(glue::glue('{functions_directory}{i}'))
  
}

#################################################################
# Import TPM values from SalmonTE
#################################################################

files = list.files(path="~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/count_tables/SalmonTE", pattern="EXPR*", full.names=TRUE, recursive=FALSE)

counter = 0
for (file in files){
  
  input = read.csv(file = file)
  if (counter == 0){TPM = input}
  else{TPM = merge(TPM, input, by = 'TE')}
  counter = counter + 1
  
}

row.names(TPM) = TPM$TE
TPM = select(TPM, -TE)
TPM = TPM[rowSums(TPM) >= 2, ]

for (i in 1:length(colnames(TPM))){
  
  if (grepl('ESC', colnames(TPM)[i])){
    
    current_name = colnames(TPM)[i]
    colnames(TPM)[i] = glue::glue('{current_name}_SRA')
    
  }
  
  if (grepl('mTEC', colnames(TPM)[i])){
    
    current_name = colnames(TPM)[i]
    colnames(TPM)[i] = glue::glue('{current_name}_ours')
    
  }
  
  else if (!(grepl('ESC', colnames(TPM)[i]))){
    
    current_name = colnames(TPM)[i]
    
    colnames(TPM)[i] = glue::glue('{current_name}_GTEx')
    
  }
  
}

#################################################################
# Batch correction
#################################################################

ID = colnames(TPM)
sampleInfo = data.frame(ID,row.names=colnames(TPM))
sampleInfo = suppressWarnings(separate(sampleInfo, col = ID, into = c('patient', 'tissue', 'batch'), sep = '_'))
sampleInfo$ID = row.names(sampleInfo)

#TPM = limma::removeBatchEffect(TPM, sampleInfo$batch)TPM

#################################################################
# PCA
#################################################################

#################################################################
# Heatmap
#################################################################

## Collapse tissue replicates

counter = 1
for (i in (unique(sampleInfo$tissue))){

  samples = subset(sampleInfo, tissue == i)$ID
  
  entry = rowMeans(TPM[ , colnames(TPM) %in% samples])
  
  if (counter == 1){
    
    output = data.frame('placeholder' = entry)
    names(output) = i
    
  }
  
  else{
    
    output[i] = entry
    
  }
  
  counter = counter + 1
  
}

averaged_TPM = output

## Rename tissues

for (name in 1:length(colnames(averaged_TPM))){
  
  if(colnames(averaged_TPM)[name] == 'Small.Intestine.Terminal.Ileum'){
    
    colnames(averaged_TPM)[name] = 'Small intestine'
    
  }
  
  
  if(colnames(averaged_TPM)[name] == 'Adrenal.Gland'){
    
    colnames(averaged_TPM)[name] = 'Adrenal gland'
    
  }
  
  if(colnames(averaged_TPM)[name] == 'Brain.Cerebellum'){
    
    colnames(averaged_TPM)[name] = 'Cerebellum'
    
  }
  
  if(colnames(averaged_TPM)[name] == 'ESC'){
    
    colnames(averaged_TPM)[name] = 'ESC'
    
  }
  
  if(colnames(averaged_TPM)[name] == 'mTEC.hi'){
    
    colnames(averaged_TPM)[name] = 'mTEC-hi'
    
  }
  
  if(colnames(averaged_TPM)[name] == 'mTEC.lo'){
    
    colnames(averaged_TPM)[name] = 'mTEC-lo'
    
  }
  
  if(colnames(averaged_TPM)[name] == 'Colon.Transverse'){
    
    colnames(averaged_TPM)[name] = "Colon"
    
  }
  
  if(colnames(averaged_TPM)[name] == 'Muscle.Skeletal'){
    
    colnames(averaged_TPM)[name] = "Skeletal muscle"
    
  }
  
  if(colnames(averaged_TPM)[name] == 'Breast.Mammary.Tissue'){
    
    colnames(averaged_TPM)[name] = "Breast"
    
  }
  
  if(colnames(averaged_TPM)[name] == 'Esophagus.Mucosa'){
    
    colnames(averaged_TPM)[name] = "Esophagus"
    
  }
  
  if(colnames(averaged_TPM)[name] == 'Adipose.Subcutaneous'){
  
    colnames(averaged_TPM)[name] = "Adipose"
    
  }
  
  if(colnames(averaged_TPM)[name] == 'Brain.Substantia.nigra'){
    
    colnames(averaged_TPM)[name] = "Substantia nigra"
    
  }
  
}

## Plot heatmap

my_heatmap = pheatmap(averaged_TPM, 
                      cluster_rows=T,
                      show_rownames=F,
                      show_colnames = T,
                      cluster_cols=T,
                      scale = 'row',
                      angle_col = 45)

save_pheatmap_png(x = my_heatmap, 
                  filename = "~/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/Plots/TPM_heatmap.png",
                  width = 3000,
                  height = 1500,
                  res = 300)

#################################################################
# Fraction of reads mapping to TEs 
#################################################################
