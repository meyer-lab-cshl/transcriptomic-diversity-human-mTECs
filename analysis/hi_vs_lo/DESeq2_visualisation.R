library(ggplot2)

##This is a set of plot templates for visualizing DESeq2 output. It assumes
##the following assignments:
##
## DESeq dataset object: dds
## Differential expression results: res
## Transformed counts: vs_dds

plotPCA(vs_dds, intgroup = 'group') + theme_bw()