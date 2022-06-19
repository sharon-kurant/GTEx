install.packages("glmnet")
library(glmnet)

# newdf = data.frame(cbind(modelX, yframe[,2]))
#colnames(newdf)[length(colnames(newdf))]="label"
#split1 = sample(c(rep(0, 0.7 * nrow(newdf)), rep(1, 0.3 * nrow(newdf))))
#train = newdf[split1==0,]
#test = newdf[split1==1,]

x <- data.matrix(modelX) # all X vars
y <- as.double(yframe[,2]) # Only Class

# Fit the LASSO model (Lasso: Alpha = 1)
#set.seed(100)
#cv.lasso <- cv.glmnet(x, y, family='multinomial', alpha=1, parallel=TRUE, standardize=TRUE, type.measure='auc')

# Results
#plot(cv.lasso)

eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}

lambdas <- 10^seq(2, -3, by = -.1)

# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(x, y, alpha = 1, standardize = TRUE, nfolds = 5)
#plot(lasso_reg, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylab = "AUC")

# Best 
#lambda_best <- lasso_reg$lambda.min 
#lambda_best


df_coef <- round(as.matrix(coef(lasso_reg, s=lasso_reg$lambda.min)), 2)
df_coef[df_coef[, 1] != 0, ]

#lasso_model <- glmnet(x, y, alpha = 1, lambda = lambda_best, standardize = TRUE)


#predictions_train <- predict(lasso_model, s = lambda_best, newx = x)
#eval_results(y, predictions_train, modelX)

cat('Min Lambda: ', lasso_reg$lambda.min, '\n 1Sd Lambda: ', lasso_reg$lambda.1se)
df_coef <- round(as.matrix(coef(lasso_reg, s=lasso_reg$lambda.min)), 2)
coef_genes = names(df_coef[which(abs(df_coef)>0),])[-1]

yframe=t(qn.brain.res.new[generankdf$genename[1:10],])

t(qn.colon.res.braincomp)[,coef_genes]
