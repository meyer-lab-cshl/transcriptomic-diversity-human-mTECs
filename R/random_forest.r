library(randomForest)
library(ROCR)
library(dplyr)
library(ggplot2)


###read in data sets and shuffle them randomly
total_int <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_internal_data_mm9.csv", sep=",", header = TRUE)
total_int <- total_int[sample(nrow(total_int)),]
write.table(total_int, file ="/home/stroemic/hiwi_16/analysis/r_random-forest/features_internal_data_mm9_for_rf.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
total_fa <- read.table("/home/stroemic/hiwi_16/analysis/r_random-forest/features_fantom5_data_mm9.csv", sep=",", header = TRUE)
total_fa <- total_fa[sample(nrow(total_fa)),]
write.table(total_fa, file ="/home/stroemic/hiwi_16/analysis/r_random-forest/features_fantom5_data_mm9_for_rf.csv", row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

#create tables only containing features 
features_int <- total_int[8:33]
#features_int <- dplyr::select(features_int, -c(Basesafter))
features_fa <- total_fa[8:33]
#features_fa <- dplyr::select(features_fa, -c(Basesafter))


folds_int <- cut(seq(1,nrow(total_int)),breaks=10,labels=FALSE)
folds_fa <- cut(seq(1,nrow(total_fa)),breaks=10,labels=FALSE)

###Retrieve general RF information
#find best mtry cutoff
bestmtry <- tuneRF(x=features_int,y=as.factor(total_int$Fantom5), ntreeTry = 100, trace = TRUE, plot = TRUE, doBest = FALSE)
#train rf without further specifications
count_no <- nrow(subset(total_int, Fantom5=='no' ))
count_big <- 2*count_no
fit_int<- randomForest(x=features_int,y=as.factor(total_int$Fantom5),
                       mtry = bestmtry[bestmtry[, 2] == min(bestmtry[, 2]), 1],
                       sampsize = c(count_no,count_big),
                       strata = as.factor(total_int$Fantom5),
                       importance = TRUE)
saveRDS(fit_int,  file="/home/stroemic/hiwi_16/analysis/r_random-forest/int_fit_s.rds")
print("saved fit_int")

#retrieve information about feature importance and plot 
importance <- data.frame(lables=row.names(fit_int$importance),Accuracy=fit_int$importance[,3],SD=fit_int$importanceSD[,3])
pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/int_feature_importance.pdf")
print(ggplot(importance, aes(x=lables, y=Accuracy)) + 
        geom_bar(stat = 'identity', fill = '#2171b5') + 
        geom_errorbar(aes(ymin=Accuracy-SD, ymax=Accuracy+SD)) +
        labs(title="Mean decrease in accuracy of RF based on internal data", x = "Features", y ="Mean Decrease Accuracy" ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
)
print(varImpPlot(fit_int, sort = TRUE, main="Importance of features in RF based on internal data"))
dev.off()

#ROc calculations
pred_int = predict(fit_int, type='prob')
pred_roc_int <- prediction(pred_int[,2], as.factor(total_int$Fantom5))
auc_int=performance(pred_roc_int, "auc")
roc_int=performance(pred_roc_int, "tpr", "fpr")
tpr_tnr_int=performance(pred_roc_int, "tpr", "tnr")



#Perform 10 fold cross validation
fits <- list()
predictions <- list()
predictions_roc <- list()

for(i in 1:10){
  print(i)
  #Segement your data by fold using the which() function
  testIndexes <- which(folds_int==i,arr.ind=TRUE)
  testData <- features_int[testIndexes, ]
  trainData <- features_int[-testIndexes, ]
  count_no <- nrow(subset(total_int[-testIndexes,], Fantom5=='no' ))
  count_big <- 2*count_no
  #Use the test and train data partitions however you desire...
  bestmtry <- tuneRF(x=trainData,y=as.factor(total_int$Fantom5[-testIndexes]), ntreeTry = 100, trace = TRUE, plot = FALSE, doBest = FALSE)
  fit <- randomForest(x=trainData,y=as.factor(total_int$Fantom5[-testIndexes]),
                      mtry = bestmtry[bestmtry[, 2] == min(bestmtry[, 2]), 1],
                      sampsize = c(count_no,count_big),
                      strata = as.factor(total_int$Fantom5[-testIndexes]),
                      importance = FALSE)
  pred <- predict(fit, newdata=testData, type='prob')
  pred_roc <- prediction(pred[,2], as.factor(total_int$Fantom5[testIndexes]))
  fits[[i]] <- fit
  predictions[[i]] <- pred
  predictions_roc[[i]] <- pred_roc
}

saveRDS(c(fits, predictions, predictions_roc), "/home/stroemic/hiwi_16/analysis/r_random-forest/int_lists_fits_predictions_predictions-roc_10x.rds")
temp1 <- readRDS("/home/stroemic/hiwi_16/analysis/r_random-forest/int_lists_fits_predictions_predictions-roc_10x.rds")
predictions_roc <- temp1[21:30]

rocs <- lapply(predictions_roc, function(x) performance(x, "tpr","fpr"))
aucs <- lapply(predictions_roc, function(x) performance(x, "auc"))
tpr_tnrs <- lapply(predictions_roc, function(x) performance(x, "tpr", "tnr"))


pdf("/home/stroemic/hiwi_16/analysis/r_random-forest/plots/int_ROC_TPR_TNR_cross_val.pdf")
plot(roc_int, main="ROC Curve for RF based on internal data",col=2,lwd=3)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
for (i in rocs){
  lines(i@x.values[[1]],i@y.values[[1]], col=4,lwd=1)
}
legend("bottomright", legend=c("Training data", "10x cross validations"),lwd = 2, col=c(2,4))

plot(tpr_tnr_int, main="TPR over TNR for RF based on internal data",col=2,lwd=3)
abline(a=1,b=-1,lwd=2,lty=2,col="gray")
for (i in tpr_tnrs){
  lines(i@x.values[[1]],i@y.values[[1]], col=4,lwd=1)
}
legend("bottomleft", legend=c("Training data", "10x cross validations"),lwd = 2, col=c(2,4))

dev.off()


######Fantom5 RF
###Retrieve general RF information
#find best mtry cutoff
bestmtry <- tuneRF(x=features_fa,y=as.factor(total_fa$Internal), ntreeTry = 100, trace = TRUE, plot = FALSE, doBest = FALSE)
#train rf without further specifications
count_no <- nrow(subset(total_fa, Internal=='no' ))
count_big <- 2*count_no
fit_fa<- randomForest(x=features_fa,y=as.factor(total_fa$Internal),
                       mtry = bestmtry[bestmtry[, 2] == min(bestmtry[, 2]), 1],
                       sampsize = c(count_no,count_big),
                       strata = as.factor(total_fa$Internal),
                       importance = TRUE)
saveRDS(fit_fa,  file="/home/stroemic/hiwi_16/analysis/r_random-forest/fa_fit_s.rds")
print("saved fit_fa")

#retrieve information about feature importance and plot 
importance <- data.frame(lables=row.names(fit_fa$importance),Accuracy=fit_fa$importance[,3],SD=fit_fa$importanceSD[,3])
pdf(file="/home/stroemic/hiwi_16/analysis/r_random-forest/plots/fa_feature_importance.pdf")
print(ggplot(importance, aes(x=lables, y=Accuracy)) + 
        geom_bar(stat = 'identity', fill = '#2171b5') + 
        geom_errorbar(aes(ymin=Accuracy-SD, ymax=Accuracy+SD)) +
        labs(title="Mean decrease in accuracy of RF based on fantom5 data", x = "Features", y ="Mean Decrease Accuracy" ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
)
print(varImpPlot(fit_fa, sort = TRUE, main="Importance of features in RF based on fantom5 data"))
dev.off()

#ROc calculations
pred_fa = predict(fit_fa, type='prob')
pred_roc_fa <- prediction(pred_fa[,2], as.factor(total_fa$Internal))
auc_fa=performance(pred_roc_fa, "auc")
roc_fa=performance(pred_roc_fa, "tpr", "fpr")
tpr_tnr_fa=performance(pred_roc_fa, "tpr", "tnr")


#Perform 10 fold cross validation
fits <- list()
predictions <- list()
predictions_roc <- list()

for(i in 1:10){
  print(i)
  #Segement your data by fold using the which() function
  testIndexes <- which(folds_fa==i,arr.ind=TRUE)
  testData <- features_fa[testIndexes, ]
  trainData <- features_fa[-testIndexes, ]
  count_no <- nrow(subset(total_fa[-testIndexes,], Internal=='no' ))
  count_big <- 2*count_no
  #Use the test and train data partitions however you desire...
  bestmtry <- tuneRF(x=trainData,y=as.factor(total_fa$Internal[-testIndexes]), ntreeTry = 100, trace = TRUE, plot = FALSE, doBest = FALSE)
  fit <- randomForest(x=trainData,y=as.factor(total_fa$Internal[-testIndexes]),
                      mtry = bestmtry[bestmtry[, 2] == min(bestmtry[, 2]), 1],
                      sampsize = c(count_no,count_big),
                      strata = as.factor(total_fa$Internal[-testIndexes]),
                      importance = FALSE)
  pred <- predict(fit, newdata=testData, type='prob')
  pred_roc <- prediction(pred[,2], as.factor(total_fa$Internal[testIndexes]))
  fits[[i]] <- fit
  predictions[[i]] <- pred
  predictions_roc[[i]] <- pred_roc
}

saveRDS(c(fits, predictions, predictions_roc), "/home/stroemic/hiwi_16/analysis/r_random-forest/fa_lists_fits_predictions_predictions-roc_10x.rds")
temp2 <- readRDS("/home/stroemic/hiwi_16/analysis/r_random-forest/fa_lists_fits_predictions_predictions-roc_10x.rds")
predictions_roc <- temp2[21:30]

rocs <- lapply(predictions_roc, function(x) performance(x, "tpr","fpr"))
aucs <- lapply(predictions_roc, function(x) performance(x, "auc"))
tpr_tnrs <- lapply(predictions_roc, function(x) performance(x, "tpr", "tnr"))


pdf("/home/stroemic/hiwi_16/analysis/r_random-forest/plots/fa_ROC_TPR_TNR_cross_val.pdf")
plot(roc_fa, main="ROC Curve for RF based on fantom5 data",col=2,lwd=3)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
for (i in rocs){
  lines(i@x.values[[1]],i@y.values[[1]], col=4,lwd=1)
}
legend("bottomright", legend=c("Training data", "10x cross validations"),lwd = 2, col=c(2,4))

plot(tpr_tnr_fa, main="TPR over TNR for RF based on fantom5 data",col=2,lwd=3)
abline(a=1,b=-1,lwd=2,lty=2,col="gray")
for (i in tpr_tnrs){
  lines(i@x.values[[1]],i@y.values[[1]], col=4,lwd=1)
}
legend("bottomleft", legend=c("Training data", "10x cross validations"),lwd = 2, col=c(2,4))

dev.off()























###########################################################################
# 
# ###determine which specifications improve TNR 
# bestmtry <- tuneRF(x=features_int_90,y=as.factor(total_int_90$Fantom5), ntreeTry = 100, trace = TRUE, plot = TRUE, doBest = FALSE)
# fit_int<- randomForest(x=features_int_90,y=as.factor(total_int_90$Fantom5),
#                         mtry = bestmtry[bestmtry[, 2] == min(bestmtry[, 2]), 1],
#                         # sampsize = c(1353,1353),
#                         # classwt = c(0.9,0.1),
#                         # strata = as.factor(total_int_90$Fantom5),
#                         importance = TRUE)
# 
# pred = predict(fit_int, newdata=features_int_10, type='prob')
# pred_roc <- prediction(pred[,2], as.factor(total_int_10[[6]]))
# auc_int=performance(pred_roc, "auc")
# pred_int=performance(pred_roc, "tpr", "fpr")
# 
# pred = predict(fit_int_c, newdata=features_int_10, type='prob')
# pred_roc <- prediction(pred[,2], as.factor(total_int_10[[6]]))
# auc_c=performance(pred_roc, "auc")
# pred_c=performance(pred_roc, "tpr", "fpr")
# 
# pred = predict(fit_int_s, newdata=features_int_10, type='prob')
# pred_roc <- prediction(pred[,2], as.factor(total_int_10[[6]]))
# auc_s=performance(pred_roc, "auc")
# pred_s=performance(pred_roc, "tpr", "fpr")
# 
# pred = predict(fit_int_all, newdata=features_int_10, type='prob')
# pred_roc <- prediction(pred[,2], as.factor(total_int_10[[6]]))
# auc_all=performance(pred_roc, "auc")
# pred_all=performance(pred_roc, "tpr", "fpr")
# 
# pdf("/home/stroemic/hiwi_16/analysis/r_random-forest/plots/ROCs.pdf")
# plot(pred_int,main="ROC Curve for Random Forest",col=2,lwd=2)
# abline(a=0,b=1,lwd=2,lty=2,col="gray")
# 
# plot(pred_c,main="ROC Curve for Random Forest with class weights",col=2,lwd=2)
# abline(a=0,b=1,lwd=2,lty=2,col="gray")
# 
# plot(pred_s,main="ROC Curve for Random Forest with sample sizes",col=2,lwd=2)
# abline(a=0,b=1,lwd=2,lty=2,col="gray")
# 
# plot(pred_all,main="ROC Curve for Random Forest with both corrections",col=2,lwd=2)
# abline(a=0,b=1,lwd=2,lty=2,col="gray")
# 
# dev.off()
# 
# 
