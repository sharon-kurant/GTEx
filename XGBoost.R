#XGBOOST
#install.packages("drat", repos="https://cran.rstudio.com")
#drat:::addRepo("dmlc"):
#  install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")
install.packages("InformationValue")#install.packages("caret")

library(caret)
library(xgboost)
library(InformationValue)
library(glmnet)

crazy_xgboost<-function(modelX, gene_index){
  # newdf -> dataML
  dataML = data.frame(cbind(modelX, yframe[,gene_index]))
  
  set.seed(100)  # For reproducibility
  # Create index for testing and training data
  inTrain <- createDataPartition(y = dataML[,length(dataML)], p = 0.8, list = FALSE)
  # subset power_plant data to training
  training <- dataML[inTrain,]
  # subset the rest to test
  testing <- dataML[-inTrain,]
  
  X_train = as.matrix(training[,-length(training)])
  y_train = training[,length(training)]
  X_test = as.matrix(testing[,-length(testing)])
  y_test = testing[,length(testing)]
  
  xgb_trcontrol = trainControl(
    method = "cv",
    number = 5,  
    allowParallel = TRUE,
    verboseIter = FALSE,
    returnData = FALSE
  )
  
  # I am specifing the same parameters with the same values as I did for Python above.
  # The hyperparameters to optimize 
  xgbGrid <- expand.grid(nrounds = c(15, 20, 25),  
                         max_depth = c(10, 15, 20),
                         colsample_bytree = seq(0.5, 0.9, length.out = 3),
                         eta = 0.1,
                         gamma=0,
                         min_child_weight = 1,
                         subsample = 1
  )
  
  set.seed(0) 
  xgb_model = train(
    X_train, y_train,  
    trControl = xgb_trcontrol,
    tuneGrid = xgbGrid,
    method = "xgbTree",
    verbose = FALSE,
    verbosity = 0
  )
  
  print(xgb_model$bestTune)
  
  predicted = predict(xgb_model, X_test)
  residuals = y_test - predicted
  MSE = mean(residuals^2)
  cat('The mean square error of the test data is ', round(MSE,3),'\n')
  
  y_test_mean = mean(y_test)
  # Calculate total sum of squares
  tss =  sum((y_test - y_test_mean)^2)
  # Calculate residual sum of squares
  rss =  sum(residuals^2)
  # Calculate R-squared
  rsq  =  1 - (rss/tss)
  cat('The R-square of the test data is ', round(rsq,3), '\n')
  
  #options(repr.plot.width=8, repr.plot.height=4)
  #my_data = as.data.frame(cbind(predicted = predicted,
  #                              observed = y_test))
  # Plot predictions vs test data
  #ggplot(my_data,aes(predicted, observed)) + geom_point(color = "darkred", alpha = 0.5) + 
  #  geom_smooth(method=lm)+ ggtitle('Linear Regression ') + ggtitle("Extreme Gradient Boosting: Prediction vs Test Data") +
  #  xlab("Predecited Power Output ") + ylab("Observed Power Output") + 
  #  theme(plot.title = element_text(color="darkgreen",size=16,hjust = 0.5),
  #        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12,hjust=.5),
  #        axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))
  
  return(list("mse" = MSE, "cor" = rsq))
}

mselist=rep(NA,10)
corrlist=rep(NA,10)

lasso_mselist=rep(NA,10)
lasso_corrlist=rep(NA,10)

xgbcorr_mselist=rep(NA,10)
xgbcorr_corrlist=rep(NA,10)


for(i in 1:10){

  mse_cor = crazy_xgboost(modelX, i)
  mselist[i]=mse_cor$mse
  corrlist[i]=mse_cor$cor
  
  
  #lasso
  lambdas <- 10^seq(2, -3, by = -.1)
  
  # Setting alpha = 1 implements lasso regression
  lasso_reg <- cv.glmnet(modelX, yframe[,i], alpha = 1, standardize = TRUE, nfolds = 5)
  #cat('Min Lambda: ', lasso_reg$lambda.min, '\n 1Sd Lambda: ', lasso_reg$lambda.1se)
  df_coef <- round(as.matrix(coef(lasso_reg, s=lasso_reg$lambda.min)), 2)
  coef_genes = names(df_coef[which(abs(df_coef)>0),])[-1]
  crazy_xgboost(modelX[,coef_genes], i)
  #x = X_train[,coef_genes]
  #xgboost 2
  mse_cor_lasso = crazy_xgboost(modelX[,coef_genes], i)
  lasso_mselist[i]=mse_cor_lasso$mse
  lasso_corrlist[i]=mse_cor_lasso$cor
  
  #xgb with correlated genes
  mse_cor_corr = crazy_xgboost(modelX[,connectgeneslists[[i]]], i)
  xgbcorr_mselist[i]=mse_cor_corr$mse
  xgbcorr_corrlist[i]=mse_cor_corr$cor
}


resultsdf1=data.frame(generankdf$genedescription[1:10])
rownames(resultsdf1)=colnames(yframe)
colnames(resultsdf1)=c("description")
resultsdf1$mse=mselist
resultsdf1$rsq=corrlist
resultsdf1$after_lasso_mse=lasso_mselist
resultsdf1$after_lasso_rsq=lasso_corrlist
resultsdf1$xgb_corr_mse=xgbcorr_mselist
resultsdf1$xgb_corr_rsq=xgbcorr_corrlist
write.csv(resultsdf1, "resultsdf_fixed.csv")


gene.f_df_fkb = as.data.frame(gene.f)
genes_fkbp = gene.f_df_fkb[gene.f_df_fkb$Name %in% coef_genes,]
write.csv(genes_fkbp, "FKBP5_lasso_genes.csv")

gene.f_df_sloc = as.data.frame(gene.f)
genes_sloc = gene.f_df_sloc[gene.f_df_sloc$Name %in% coef_genes,]
write.csv(genes_sloc, "SLCO4A1_lasso_genes.csv")


which(pheno.f$Name %in% coef_genes)
gene.f[,1]
gene.f[gene.f$]

xcorr = train[,connectgeneslists[[i]]]
mse_cor = run_xgboost(train[,connectgeneslists[[i]]],test[,connectgeneslists[[i]]],i)
xgbcorr_mselist[i]=mse_cor$mse
xgbcorr_corrlist[i]=mse_cor$cor


















set.seed(100)
max_depth = 4
num_rounds = 8

run_xgboost <- function(train, test, i){                   data               label                                                                                  "binary:logistic"
  bstSparse <- xgboost(data = as.matrix(train[,-ncol(train)]), label = as.numeric(train[,ncol(train)]), max.depth = max_depth, eta = 1.5, nthread = 2, nrounds = num_rounds, objective = "reg:squarederror")
  pred <- predict(bstSparse, as.matrix(test[,-ncol(test)]))
  
  pred
  #                      test labels
  err <- mean(pred - test[,ncol(test)])^2
  #corrlist[i]=cor(pred,test[,ncol(test)])
  
  print(paste("mse=", err))
  #mselist[i]=err
  
  return(list("mse" = err, "cor" = cor(pred,test[,ncol(test)])))
}



mselist=rep(NA,10)
corrlist=rep(NA,10)

lasso_mselist=rep(NA,10)
lasso_corrlist=rep(NA,10)

xgbcorr_mselist=rep(NA,10)
xgbcorr_corrlist=rep(NA,10)

split1 <- sample(c(rep(0, 0.7 * nrow(newdf)), rep(1, 0.3 * nrow(newdf))))

for(i in 1:10){
  newdf = data.frame(cbind(modelX, yframe[,i]))
  colnames(newdf)[length(colnames(newdf))]="label"
  train = newdf[split1==0,]
  test = newdf[split1==1,]
  
  mse_cor = run_xgboost(train,test,i)
  mselist[i]=mse_cor$mse
  corrlist[i]=mse_cor$cor
  
  #lasso
  lambdas <- 10^seq(2, -3, by = -.1)
  
  # Setting alpha = 1 implements lasso regression
  lasso_reg <- cv.glmnet(modelX, yframe[,i], alpha = 1, standardize = TRUE, nfolds = 5)
  #cat('Min Lambda: ', lasso_reg$lambda.min, '\n 1Sd Lambda: ', lasso_reg$lambda.1se)
  df_coef <- round(as.matrix(coef(lasso_reg, s=lasso_reg$lambda.min)), 2)
  coef_genes = names(df_coef[which(abs(df_coef)>0),])[-1]
  x = train[,coef_genes]
  #xgboost 2
  mse_cor = run_xgboost(train[,coef_genes],test[,coef_genes],i)
  lasso_mselist[i]=mse_cor$mse
  lasso_corrlist[i]=mse_cor$cor
  
  #xgb with correlated genes
  xcorr = train[,connectgeneslists[[i]]]
  mse_cor = run_xgboost(train[,connectgeneslists[[i]]],test[,connectgeneslists[[i]]],i)
  xgbcorr_mselist[i]=mse_cor$mse
  xgbcorr_corrlist[i]=mse_cor$cor
}
resultsdf=data.frame(generankdf$genedescription[1:10])
rownames(resultsdf)=colnames(yframe)
colnames(resultsdf)=c("description")
resultsdf$mse=mselist
resultsdf$rsq=corrlist*corrlist
resultsdf$after_lasso_mse=lasso_mselist
resultsdf$after_lasso_rsq=lasso_corrlist*lasso_corrlist
resultsdf$xgb_corr_mse=xgbcorr_mselist
resultsdf$xgb_corr_rsq=xgbcorr_corrlist*xgbcorr_corrlist
