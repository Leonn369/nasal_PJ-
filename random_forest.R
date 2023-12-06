###random forest analysis was performed to explore the discriminatory potential of sex in the nasal microbiome

library(randomForest)
library(caret)
library(tidyverse)
library(pROC)
library(e1071)
library(ggplot2)

##meta data input
meta <- read.table("./nasl_meta.tsv",sep = "\t", header = T)
meta.class <- select(meta,"Sex")

n <- c("bacteria","fungi","total")
#pred <- list()
for (i in 1:3) {
  df <- read.table(paste0(n[i], "_occ0.tsv"), sep = "\t", header = T)
  df <- df[,which(colMeans(df > 0) > 0.1)]
  df <- cbind(df, meta)
  df <- df %>% df(Sex)
  df$Sex <- factor(df$Sex)

  feature <- df[, -ncol(df)]
  label <- df[, ncol(df)]

  mtry <- floor(sqrt(ncol(feature)))
  res <- get_RF(feature = feature, label = label, mtry = mtry)
  
  assign(res, n[1])
}

##for example for bacteria result
auc <- get_AUC(bacteria, n = 100)

#chose the mean auc model in 10 fold 10 times random forest to get ROC plot
mean <- mean(auc)
data.pred1 <- data.pred[[9]] ##the mean AUC model in result list
ran_roc <- roc(data.pred1$true, data.pred1$prob,percent = TRUE,plot=TRUE,ci=TRUE)
rocobj <- plot.roc(data.pred1$true, data.pred1$prob,
                   main="Confidence intervals", 
                   percent=TRUE,ci=TRUE
) 
ciobj <- ci.se(rocobj,
               specificities=seq(0, 100, 5)
)
plot(ciobj, type="shape", col="#1c61b6AA")



get_RF <- function(feature,label,mtry){
  
  data2 <- list()
  auc_values <- c()
  
  for (i in 1:10) {
    
    folds <- createFolds(label, k=10 ,list = TRUE, returnTrain = TRUE)

    for (j in 1:10) {

      train_features <- feature[folds[[j]], ]
      train_labels <- label[folds[[j]]]
      test_features <- feature[-folds[[j]], ]
      test_labels <- label[-folds[[j]]]

      rf_model <- randomForest(x = train_features, y = train_labels,ntree = 500, mtry = mtry)
      rf_pred <- predict(rf_model, test_features,type = "prob")
      auc_value <- roc(test_labels, rf_pred[,2])$auc[1]
      data2[[10*(i-1)+j]] <- data.frame(true=test_labels, prob=rf_pred[,2],auc=auc_value)
    }
  }
  return(data2)
}

get_AUC <- function(data,n){
  auc_values <- c()
  for(i in 1:n){
    data1 <- data[[i]]
    auc1 <- unique(data1$auc)
    auc_values <- c(auc_values,auc1)
  }
  return(auc_values)
}


