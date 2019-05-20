library(plyr)
library(data.table)
library(ggplot2)
library(LSD)

setwd("/home/stroemic/hiwi_16/analysis/lars_perl/mouse/")

files <- list.files(pattern="*trimmed.isoforms.csv")
fantoms <- list.files(pattern="*ctss.isoforms.csv")
datalist=list()
c <- 1
for(i in files){
  dat <- read.table(i, header = TRUE, sep = ",")
  dat$iteration <- c
  datalist[[c]] <- dat
  c <- c+1
}

dat_internal = do.call(rbind, datalist)

datalist=list()
c <- 1
for(i in fantoms){
  dat <- read.table(i, header = TRUE, sep = ",")
  dat$iteration <- c
  datalist[[c]] <- dat
  c <- c+1
}

dat_fantoms = do.call(rbind, datalist)


dat_agg_internal <- setDT(dat_internal)[, lapply(.SD, sum), by=.(geneID, Position), .SDcols=c("BarcodeCount","ReadCount","PosFromAnno","Class")]
setDF(dat_agg_internal)
dat_agg_fantoms <- setDT(dat_fantoms)[, lapply(.SD, sum), by=.(geneID, Position), .SDcols=c("BarcodeCount","ReadCount","PosFromAnno","Class")]
setDF(dat_agg_fantoms)

write.table(dat_agg_internal, file ="/home/stroemic/hiwi_16/analysis/lars_perl/mouse/all_internal.isoforms.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(dat_agg_fantoms, file ="/home/stroemic/hiwi_16/analysis/lars_perl/mouse/all_fantoms.isoforms.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

total <- merge(dat_agg_internal, dat_agg_fantoms, by=c("geneID","Position"), all=TRUE)
total <- subset(total, select = -c(ReadCount.x, PosFromAnno.x, Class.x, ReadCount.y, PosFromAnno.y, Class.y))
total <- rename(total,c("BarcodeCount.x"="CountInternal","BarcodeCount.y"="CountFantom"))

pdf("/home/stroemic/hiwi_16/analysis/lars_perl/mouse/plots/Counts_per_position.pdf")
heatscatter(total$CountInternal, total$CountFantom, log = "xy", 
            main = "Counts per position",
            xlab = "Count internal data",
            ylab = "Count Fantom5 data")

ggplot(total, aes(x= CountInternal, y= CountFantom)) + 
  geom_point() + 
  labs(title= "Counts per position", x= "Count internal data", y= "Count Fantom5 data") +
  scale_x_log10() + 
  scale_y_log10()
dev.off()
# 
# pdf("/home/stroemic/hiwi_16/analysis/lars_perl/mouse/plots/Counts_Histogramm.pdf")
# hist(total$CountInternal, breaks=((min(total$CountInternal)-1):(max(total$CountInternal)+1)), 
#      main="Internal data", cex.lab= 1, xlab="Counts per position", ylab="Frequency" ,
#      col="#253494", xlim=c(0,max(total$CountInternal, na.rm = TRUE)))
# hist(total$CountFantom, breaks=((min(total$CountFantom)-1):(max(total$CountFantom)+1)), 
#      main="Fantom5 data", cex.lab= 1, xlab="Counts per position", ylab="Frequency" ,
#      col="#253494", xlim=c(0,max(total$CountFantom, na.rm = TRUE)))
# total_sub <- total[total$CountInternal > 100 & total$CountFantom > 100, ]
# hist(total_sub$CountInternal, breaks=((min(total_sub$CountInternal)-1):(max(total_sub$CountInternal)+1)), 
#      main="Internal data >100 Counts per Pos", cex.lab= 1, xlab="Counts per position", ylab="Frequency" ,
#      col="#253494", xlim=c(0,max(total_sub$CountInternal, na.rm = TRUE)))
# hist(total_sub$CountFantom, breaks=((min(total_sub$CountFantom)-1):(max(total_sub$CountFantom)+1)), 
#      main="Fantom5 data >100 Counts per Pos", cex.lab= 1, xlab="Counts per position", ylab="Frequency" ,
#      col="#253494", xlim=c(0,max(total_sub$CountFantom, na.rm = TRUE)))
# dev.off()