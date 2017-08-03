library(randomForest)
library(dplyr)

setwd("/home/stroemic/hiwi_16/analysis/lars_perl/all_63/features/")

fit_int <- readRDS("/home/stroemic/hiwi_16/analysis/r_random-forest/int_fit_s.rds")
fit_fa <- readRDS("/home/stroemic/hiwi_16/analysis/r_random-forest/fa_fit_s.rds")


files <- list.files(pattern = "features.csv$")

for (i in files[5:6]) {
  dat <- read.table(i, sep = ",", header=TRUE)
  features <- dat[7:32]
  if (grepl("_", i)) {
    name= paste((head(unlist(strsplit(i, "[_]" )), n=-1)), collapse='_')
    print(name)
    print("using fit_int")
    dat$prediction <- predict(fit_int, newdata = features)
  } else {
    name= paste((head(unlist(strsplit(i, "[.]" )), n=-2)), collapse='.')
    print(name)
    print("using fit_fa")
    dat$prediction <- predict(fit_fa, newdata = features)
  }
  dat_filtered <- dat[dat$prediction == 'yes', ]

  
  write.table(dat_filtered[,c(1:6)], file =paste("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/",name, ".filtered.csv", sep=""), row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
  write.table(dat_filtered, file =paste("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/",name, ".features-filtered.csv", sep=""), row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
  
  print(paste("/home/stroemic/hiwi_16/analysis/r_random-forest/filtered_data/",name, ".filtered.csv", sep=""))
  
}
