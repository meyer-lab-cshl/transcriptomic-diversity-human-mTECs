#################################################################
# PCA
#################################################################

PCA = plotPCA(vs_dds, intgroup = 'group') + 
  ggtitle('mTEC-hi vs mTEC-lo', 'TE expression PCA') +
  geom_text(aes(label = colnames(vs_dds)), nudge_x = 0.5, nudge_y = 0.2)

PCA + theme_bw() + theme(plot.title = element_text(face = 'bold', size = 20),
                         plot.subtitle = element_text(size = 14),
                         axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 14),
                         axis.title = element_text(size = 14),
                         axis.line = element_line(size = 0.8),
                         panel.border = element_blank())

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo/Plots/hi_vs_lo_TEs_PCA.png", 
       width = 20, height = 15, units = "cm")