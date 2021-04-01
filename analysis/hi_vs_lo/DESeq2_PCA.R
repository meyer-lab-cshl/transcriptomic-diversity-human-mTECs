library(ggplot2)
library(DESeq2)
library(limma)

#################################################################
# PCA
#################################################################

PCA = plotPCA(vs_dds_transcripts_TE, intgroup = 'tissue') + 
  ggtitle('TE expression PCA', 'TEtranscripts') 

PCA + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                         plot.subtitle = element_text(size = 14),
                         axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),
                         axis.title = element_text(size = 14),
                         axis.line = element_line(size = 0.8),
                         panel.border = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/21-04-01/PCA.png", 
       width = 20, height = 20, units = "cm")