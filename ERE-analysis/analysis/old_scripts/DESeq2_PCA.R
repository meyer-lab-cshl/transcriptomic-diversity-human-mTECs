library(ggplot2)
library(DESeq2)

#################################################################
# PCA
#################################################################

pcaData = plotPCA(vs_dds_transcripts_TE, intgroup='tissue', returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

PCA = ggplot(pcaData, aes(PC1, PC2, fill = tissue)) + 
  geom_point(size=3, shape = 21, stroke = 0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  labs(fill= "Tissue") +
  scale_fill_manual(values = c('#4c72b0ff', '#dd8452ff')) 

PCA + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                         plot.subtitle = element_text(size = 14),
                         axis.text.x = element_text(size = 11),
                         axis.text.y = element_text(size = 11),
                         axis.title.x = element_text(size = 13, margin = margin(t = 8)),
                         axis.title.y = element_text(size = 13),
                         axis.line = element_line(size = 0.8),
                         panel.border = element_blank(),
                         legend.text = element_text(size = 10),
                         legend.title = element_text(size = 13))

ggsave("/Users/mpeacey/TE_thymus/analysis/Plots/Presentation/transcripts_PCA.png", 
       width = 7, height = 6, units = "in")