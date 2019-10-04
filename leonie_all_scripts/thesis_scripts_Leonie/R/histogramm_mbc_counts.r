library(gtools)
library(dplyr)

setwd("/home/stroemic/hiwi_16/analysis/STAR_alignment_hg19/5Seq_processed_IVT/")

folders <- mixedsort(list.files(pattern ="^C.*"))
table_files=grep('collapsing_1/.*countsPerMol',list.files(path = folders, all.files = TRUE, full.names = TRUE, recursive = TRUE),value = TRUE)

pdf("/home/stroemic/hiwi_16/documentation/histogramms_mbc_counts.pdf")
for (i in table_files){
  x=sapply(strsplit(i, "/"), '[',3)
  bed=read.table(paste(getwd(),"/",i,sep=""), header = TRUE)
  bed= filter(bed, count >= 3 & count <= 5000)
  histogramm=hist(bed$count, breaks=(-0.5:(max(bed$count)+0.5)), main=x, cex.lab= 1, xlab="Counts per mbc", ylab="Frequency",col="#33a02c", xlim=c(0,max(bed$count, na.rm = TRUE)))
}
dev.off()
