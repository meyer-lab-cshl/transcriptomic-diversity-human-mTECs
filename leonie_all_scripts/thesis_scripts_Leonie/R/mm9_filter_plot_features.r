library(randomForest)
library(ROCR)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)

###read in random forest classifier
fit_int <- readRDS("/home/stroemic/hiwi_16/analysis/r_random-forest/int_fit_s.rds")
fit_fa <- readRDS("/home/stroemic/hiwi_16/analysis/r_random-forest/fa_fit_s.rds")

####filter mouse data
ms_int <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_internal_data_mm9.csv", sep = ",", header = TRUE)
ms_fa <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_fantom5_data_mm9.csv", sep = ",", header = TRUE)

ms_int_features <- ms_int[8:33]
ms_fa_features <- ms_fa[8:33]
##predict occurence of position when found with other method
ms_int$prediction <- predict(fit_int, newdata=ms_int_features)
ms_fa$prediction <- predict(fit_fa, newdata =ms_fa_features)
##filter out all positions which wouldn't have been found with both methods
ms_int_filtered <- ms_int[ms_int$prediction == 'yes',]
ms_int_filtered <- dplyr::rename(ms_int_filtered, BarcodeCount=BarcodeCount.x)
ms_fa_filtered <- ms_fa[ms_fa$prediction == 'yes', ]
ms_fa_filtered <- dplyr::rename(ms_fa_filtered, BarcodeCount=BarcodeCount.x)

write.table(ms_int_filtered, file ="/home/stroemic/hiwi_16/analysis/r_random-forest/features_internal_data_mm9.filtered.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(ms_fa_filtered, file ="/home/stroemic/hiwi_16/analysis/r_random-forest/features_fantom5_data_mm9.filtered.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(ms_int_filtered[,c(1:5,7)], file ="/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all_internal_mm9_ESC.filtered.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(ms_fa_filtered[,c(1:5,7)], file ="/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/all.fantom5.mm9.ESC.filtered.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")




###plot new plots with filtered data
##calculate dataframes only containing internal, fantom5 and both positions
#internal
ms_int_only <- merge(ms_int_filtered, ms_fa_filtered, by=c("geneID","chrom","strand", "start", "end"), all.x=TRUE)
ms_int_only <- ms_int_only[, -grep(".y$", colnames(ms_int_only))]
ms_int_only[is.na(ms_int_only$Internal),]$Fantom5<- 'no'
ms_int_only[!is.na(ms_int_only$Internal),]$Fantom5 <- 'yes'
ms_int_only <- ms_int_only[ms_int_only$Fantom5 == 'no',]
ms_int_only <- select(ms_int_only, -c(Fantom5, prediction.x, Internal))
colnames(ms_int_only) <- sub(".x$", "", colnames(ms_int_only))
ms_int_only$source <- "Internal_only"
#fantom5
ms_fa_only <- merge(ms_fa_filtered, ms_int_filtered, by=c("geneID","chrom","strand", "start", "end"), all.x=TRUE)
ms_fa_only <- ms_fa_only[, -grep(".y$", colnames(ms_fa_only))]
ms_fa_only[is.na(ms_fa_only$Fantom5),]$Internal<- 'no'
ms_fa_only[!is.na(ms_fa_only$Fantom5),]$Internal <- 'yes'
ms_fa_only <- ms_fa_only[ms_fa_only$Internal == 'no',]
ms_fa_only <- select(ms_fa_only, -c(Fantom5, prediction.x, Internal))
colnames(ms_fa_only) <- sub(".x$", "", colnames(ms_fa_only))
ms_fa_only$source <- "Fantom5_only"
#both
total_both <- merge(ms_int_filtered, ms_fa_filtered, by=c("geneID","chrom","strand", "start", "end", "prediction") )
total_both <- select(total_both, -c(Internal, Fantom5, prediction))
total_both_int <- total_both[, -grep(".y$", colnames(total_both))]
colnames(total_both_int) <- sub(".x$", "", colnames(total_both_int))
total_both_int$source <- "Both"
total_both_fa <- total_both[, -grep(".x$", colnames(total_both))]
colnames(total_both_fa) <- sub(".y$", "", colnames(total_both_fa))
total_both_fa$source <- "Both_Fantom5"

massive <- rbind(ms_int_only, ms_fa_only,  total_both_int)


###Plots features
#Freq of bases
pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/filtered_features_Freq_bases.pdf")
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
pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/filtered_features_absolut_positions.pdf")
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

pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/filtered_features_absolute_positions_counts_proportions.pdf")
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
