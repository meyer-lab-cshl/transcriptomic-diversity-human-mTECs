library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)

total_int <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_internal_data_mm9.csv", sep=",", header = TRUE)
total_fa <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_fantom5_data_mm9.csv", sep=",", header = TRUE)

total_int_only <- total_int[total_int$Fantom5 == 'no',]
total_int_only <- total_int_only[,-6]
total_int_only$source <- "Internal_only"
total_fa_only <- total_fa[total_fa$Internal == 'no',]
total_fa_only <- total_fa_only[,-6]
total_fa_only$source <- "Fantom5_only"

total_both <- merge(total_int, total_fa, by=c("geneID","chrom","strand", "start", "end") )
total_both <- total_both[,-c(6,34)]
total_both_int <- total_both[, -grep(".y$", colnames(total_both))]
colnames(total_both_int) <- sub(".x$", "", colnames(total_both_int))
total_both_int$source <- "Both"

massive <- rbind(total_int_only, total_fa_only,  total_both_int)
###Plots features
#Freq of bases
pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/features_Freq_bases.pdf")
my.x.labs <- list("Both","Fantom5 only","Internal only")
##plot Freq10 and 50 seperately for label names
for (i in grep("Freq10",colnames(massive)) ) {
  print(colnames(massive)[i])
  print(ggplot(massive, aes(x=source, y=massive[[i]], fill=source)) + 
          theme_bw() +
          geom_boxplot() + 
          labs(title=paste("Feature \"",colnames(massive)[i],"\" in different datasets", sep = ""), x ="", y =paste("Counts of base ",tail(strsplit(colnames(massive)[i],"")[[1]], n=1)," in 10bp window",sep="") ) +
          scale_fill_manual(values = c('#2171b5','#2171b5','#2171b5','#2171b5'), guide=FALSE) +
          scale_x_discrete(labels = my.x.labs) )
}
for (i in grep("Freq50",colnames(massive)) ) {
  print(colnames(massive)[i])
  print(ggplot(massive, aes(x=source, y=massive[[i]], fill=source)) + 
          theme_bw() +
          geom_boxplot() + 
          labs(title=paste("Feature \"",colnames(massive)[i],"\" in different datasets", sep = ""), x ="", y =paste("Counts of base ",tail(strsplit(colnames(massive)[i],"")[[1]], n=1)," in 50bp window",sep="")) +
          scale_fill_manual(values = c('#2171b5','#2171b5','#2171b5','#2171b5'), guide=FALSE) +
          scale_x_discrete(labels = my.x.labs) )
}
dev.off()
#absolut positions
pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/features_absolut_positions.pdf")
for (i in (10:23)) {
  print(colnames(massive)[i])
  print(ggplot(massive, aes(x=massive$source)) + 
          theme_bw() +
          geom_bar(aes(fill=massive[[i]]), position='fill') +
          labs(title=paste("Feature \"",colnames(massive)[i],"\" in different datasets", sep = ""), x ="", y =colnames(massive)[i] ) +
          scale_fill_manual(values = c('#2171b5','#084594')) +
          scale_x_discrete(labels = my.x.labs) +
          guides(fill=guide_legend(title=NULL))
        )
  
}
dev.off()

#positions stacked
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))

positions_tss_tts <- massive[,c(10:13,17,19,18,20,33)]
positions_others <- massive[,c(14:16,21:23,33)]
positions_all <- massive[,c(10:23,33)]
m_pos_tss_tts <- melt(positions_tss_tts, id.vars = "source")
m_pos_others <- melt(positions_others, id.vars = "source")
m_pos_all <- melt(positions_all, id.vars = "source")
###aggregate false and true values to counts
m_pos_tss_tts <- aggregate(m_pos_tss_tts$value, by=list(m_pos_tss_tts$source, m_pos_tss_tts$variable), FUN=sum )
m_pos_others <- aggregate(m_pos_others$value, by=list(m_pos_others$source, m_pos_others$variable), FUN=sum )
m_pos_all <- aggregate(m_pos_all$value, by=list(m_pos_all$source, m_pos_all$variable), FUN=sum )

pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/features_absolute_positions_counts_proportions.pdf")
my.x.labs <- list("Both","Fantom5 only","Internal only")
my.labs <- list("TSS \u00B1 10 genes","TSS \u00B1 10 transcripts","TSS \u00B1 100 genes","TSS \u00B1 100 transcripts","First exon",
                "Other exon","Intron","TTS \u00B1 100 genes","TTS + 200 genes", "TTS \u00B1 100 transcripts","TTS + 200 transcripts",
                "1kb downstream","1kb upstream","Antisense")
print(ggplot(m_pos_tss_tts, aes(x=Group.1,y=x,fill=Group.2)) +
        theme_bw() +
        geom_bar(stat = "identity") +
        labs(title="Number of positions belonging to each feature in different datasets", x ="", y ="Count") +
        #scale_fill_manual(values = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594')) +
        scale_fill_brewer(palette="Blues") +
        scale_x_discrete(labels = my.x.labs) +
        guides(fill=guide_legend(title="Feature"))
      )
print(ggplot(m_pos_others, aes(x=Group.1,y=x,fill=Group.2)) +
        theme_bw() +
        geom_bar(stat = "identity") +
        labs(title="Number of positions belonging to each feature in different datasets", x ="", y ="Count") +
        #scale_fill_manual(values = c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')) +
        scale_fill_brewer(palette="Blues") +
        scale_x_discrete(labels = my.x.labs) +
        guides(fill=guide_legend(title="Feature"))
      )
print(ggplot(m_pos_all, aes(x=Group.1,y=x,fill=Group.2)) + 
        theme_bw() +
        geom_bar(position="fill",stat = "identity") +
        labs(title="Proportion of positions belonging to each feature in different datasets", x ="", y ="Proportion") +
        #scale_fill_manual(values = c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')) +
        scale_fill_discrete(labels = my.labs) +
        scale_x_discrete(labels = my.x.labs) +
        guides(fill=guide_legend(title="Feature"))
)
dev.off()


# ###
# pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/features_relative_counts_positions_.pdf")
# for (i in grep("rel",colnames(massive)) ) {
#   print(colnames(massive)[i])
#   print(ggplot(massive, aes(x=source, y=massive[[i]])) + 
#           geom_boxplot() + #outlier.shape = NA
#           #scale_y_continuous(limits = quantile(massive[[i]], c(0.1, 0.9))) +
#           labs(title=paste("Feature \"",colnames(massive)[i],"\" in different datasets", sep = ""), x ="", y =colnames(massive)[i] ))
# }
# dev.off()
