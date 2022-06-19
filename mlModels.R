#XGBOOST
install.packages("drat", repos="https://cran.rstudio.com")
drat:::addRepo("dmlc"):
install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")
install.packages("xgboost")
require(xgboost)

# 
# train = as.matrix(ml_genes_df_p6[c(0:50),])
# test <- as.matrix(ml_genes_df_p6[-c(0:50),])
set.seed(1111)
max_depth = 4
num_rounds = 8
# #                             data               label                                                                                  "binary:logistic"
# bstSparse <- xgboost(data = train[,-ncol(train)], label = train[,ncol(train)], max.depth = max_depth, eta = 1.5, nthread = 2, nrounds = num_rounds, objective = "reg:squarederror")
# pred <- predict(bstSparse, test[,-ncol(test)])
# 
# 
# pred
# #                      test labels
# err <- mean(pred - test[,ncol(test)])^2
# 
# print(paste("mse=", err))
run_xgboost <- function(train, test, i){
  #newdf = data.frame(cbind(modelX, yframe[,i]))
  #colnames(newdf)[length(colnames(newdf))]="label"
  #split1 <- sample(c(rep(0, 0.7 * nrow(newdf)), rep(1, 0.3 * nrow(newdf))))
  #train = newdf[split1==0,]
  #test=newdf[split1==1,]
  #                             data               label                                                                                  "binary:logistic"
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
