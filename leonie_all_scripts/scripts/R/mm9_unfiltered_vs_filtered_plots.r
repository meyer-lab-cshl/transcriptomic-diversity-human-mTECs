library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(LSD)

ms_int <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_internal_data_mm9.csv", sep = ",", header = TRUE)
ms_fa <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_fantom5_data_mm9.csv", sep = ",", header = TRUE)

ms_int_filt <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_internal_data_mm9.filtered.csv", sep=",", header=TRUE)
ms_fa_filt <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_fantom5_data_mm9.filtered.csv", sep=",", header=TRUE)


###plot new plots with filtered data
##calculate dataframes containing internal and  fantom5 sites before and after filtering
#internal
ms_int <- select(ms_int, -c(Fantom5))
ms_int$source <- "Internal"
#fantom5
ms_fa <- select(ms_fa, -c(Internal))
ms_fa$source <- "Fantom5"
#filtered internal
ms_int_filt <- select(ms_int_filt, -c(Fantom5, prediction))
ms_int_filt <- dplyr::rename(ms_int_filt, BarcodeCount.x = BarcodeCount)
ms_int_filt$source <- "Filtered_internal"
#filtered fantom5
ms_fa_filt <- select(ms_fa_filt, -c(Internal, prediction))
ms_fa_filt <- dplyr::rename(ms_fa_filt, BarcodeCount.x = BarcodeCount)
ms_fa_filt$source <- "Filtered_fantom5"

massive <- rbind(ms_int, ms_fa, ms_int_filt, ms_fa_filt)
massive$source <- factor(massive$source, levels = c("Fantom5", "Internal", "Filtered_fantom5", "Filtered_internal"))

####Scatterplots
total <- merge(ms_int_filt[c(1,4,6)], ms_fa_filt[c(1,4,6)], by=c("geneID","start"), all=TRUE)
pdf("/home/stroemic/hiwi_16/analysis/r_random-forest/plots/Counts_per_position_filtered.pdf")
heatscatter(total$BarcodeCount.x.x, total$BarcodeCount.x.y, log = "xy", 
            main = "Counts per position",
            xlab = "Count internal data",
            ylab = "Count Fantom5 data")

ggplot(total, aes(x= BarcodeCount.x.x, y= BarcodeCount.x.y)) + 
  geom_point() + 
  labs(title= "Counts per position", x= "Count internal data", y= "Count Fantom5 data") +
  scale_x_log10() + 
  scale_y_log10()
dev.off()


###Plots features
#Freq of bases
pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/unfiltered_vs_filtered_features_Freq_bases.pdf")
for (i in grep("Freq10",colnames(massive)) ) {
  print(colnames(massive)[i])
  print(ggplot(massive, aes(x=source, y=massive[[i]], fill=source)) + 
          geom_boxplot() + 
          labs(title=paste("Feature \"",colnames(massive)[i],"\" in different datasets", sep = ""), x ="", y =paste("Counts of base ",tail(strsplit(colnames(massive)[i],"")[[1]], n=1)," in 10bp window",sep="") ) +
          scale_fill_manual(values = c('#2171b5','#2171b5','#2171b5','#2171b5'), guide=FALSE))
}
for (i in grep("Freq50",colnames(massive)) ) {
  print(colnames(massive)[i])
  print(ggplot(massive, aes(x=source, y=massive[[i]], fill=source)) + 
          geom_boxplot() + 
          labs(title=paste("Feature \"",colnames(massive)[i],"\" in different datasets", sep = ""), x ="", y =paste("Counts of base ",tail(strsplit(colnames(massive)[i],"")[[1]], n=1)," in 50bp window",sep="")) +
          scale_fill_manual(values = c('#2171b5','#2171b5','#2171b5','#2171b5'), guide=FALSE))
}
dev.off()
#absolut positions
pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/unfiltered_vs_filtered_features_absolut_positions.pdf")
for (i in (10:23)) {
  print(colnames(massive)[i])
  print(ggplot(massive, aes(x=massive$source)) + 
          geom_bar(aes(fill=massive[[i]]), position='fill') +
          labs(title=paste("Feature \"",colnames(massive)[i],"\" in different datasets", sep = ""), x ="", y =colnames(massive)[i] ) +
          scale_fill_manual(values = c('#2171b5','#084594')) +
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
m_pos_tss_tts <- aggregate(m_pos_tss_tts$value, by=list(m_pos_tss_tts$source, m_pos_tss_tts$variable), FUN=sum )
m_pos_others <- aggregate(m_pos_others$value, by=list(m_pos_others$source, m_pos_others$variable), FUN=sum )
m_pos_all <- aggregate(m_pos_all$value, by=list(m_pos_all$source, m_pos_all$variable), FUN=sum )

pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/unfiltered_vs_filtered_features_absolute_positions_counts_proportions.pdf")
print(ggplot(m_pos_tss_tts, aes(x=Group.1,y=x,fill=Group.2)) +
        geom_bar(stat = "identity") +
        labs(title="Number of positions belonging to each feature in different datasets", x ="", y ="Count") +
        #scale_fill_manual(values = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594')) +
        scale_fill_brewer(palette="Blues") +
        guides(fill=guide_legend(title="Feature"))
)
print(ggplot(m_pos_others, aes(x=Group.1,y=x,fill=Group.2)) +
        geom_bar(stat = "identity") +
        labs(title="Number of positions belonging to each feature in different datasets", x ="", y ="Count") +
        #scale_fill_manual(values = c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')) +
        scale_fill_brewer(palette="Blues") +
        guides(fill=guide_legend(title="Feature"))
)
print(ggplot(m_pos_all, aes(x=Group.1,y=x,fill=Group.2)) + 
        geom_bar(position="fill",stat = "identity") +
        labs(title="Proportion of positions belonging to each feature in different datasets", x ="", y ="Proportion") +
        #scale_fill_manual(values = c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')) +
        scale_color_manual(values=mycolors) +
        guides(fill=guide_legend(title="Feature"))
)
dev.off()
