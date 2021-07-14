library(ggvenn)

input = list(local = sigLocus_local,
             transcripts = sigLocus_transcripts)

ggvenn(input)

ggsave("/Users/mpeacey/TE_thymus/analysis/hi_vs_lo_local/Plots/local_vs_transcripts_venn.png", 
       width = 20, height = 15, units = "cm")