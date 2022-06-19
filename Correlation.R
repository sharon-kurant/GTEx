# genes = 5000
# brain_heatmap =matrix(as.numeric(t(qn.brain.res[c(0:genes),c(0:62)])),ncol =genes,dimnames=dimnames(t(qn.brain.res[c(0:genes),c(0:62)])))
# colon_heatmap =matrix(as.numeric(t(qn.colon.res[c(0:genes),])),ncol =genes,dimnames=dimnames(t(qn.colon.res[c(0:genes),])))
# 
# #colon_brain_heatmap=mapply(cor, as.data.frame(t(qn.colon.res)), as.data.frame(t(qn.brain.res[,c(0:60)])))
# colon_brain_heatmap =cor(colon_heatmap, brain_heatmap)
# #heatmap(colon_brain_heatmap)
# 
# colon_brain_cor_df = as.data.frame(colon_brain_heatmap)
# 
# colon_brain_cor_df[colon_brain_cor_df>0.6]
# x = which(colon_brain_cor_df>0.6, arr.ind = TRUE)
# 
# rownames(x)
# 
# colon_brain_heatmap[colon_brain_heatmap > 0.51]
# 
# dimnames(colon_brain_cor_df[colon_brain_cor_df > 0.51])
# 
# which.max(colon_brain_cor_df)
# colon_brain_cor_df[19,9]
# 
# colnames(colon_brain_heatmap)[apply(colon_brain_heatmap, 1, function (x) which(x==max(x[x>0.5])))]

###################### Tal & PELE Version ################# 

n_genes= 5000
brain_heatmap =matrix(as.numeric(t(qn.brain.res[c(0:n_genes),c(0:60)])),ncol =n_genes,dimnames=dimnames(t(qn.brain.res[c(0:n_genes),c(0:60)])))
colon_heatmap =matrix(as.numeric(t(qn.colon.res[c(0:n_genes),])),ncol =n_genes,dimnames=dimnames(t(qn.colon.res[c(0:n_genes),])))

#colon_brain_heatmap=mapply(cor, as.data.frame(t(qn.colon.res)), as.data.frame(t(qn.brain.res[,c(0:60)])))
colon_brain_heatmap =cor(colon_heatmap, brain_heatmap)

#heatmap(colon_brain_heatmap)
colon_brain_cor_df = as.data.frame(colon_brain_heatmap)

#which(colon_brain_cor_df > 0.45)

#colon_brain_cor_df[8,20]

#colnames(colon_brain_heatmap)[apply(colon_brain_heatmap, 1, function (x) which(x==max(x[x>0.5])))]

cols_genes_sum = colSums(abs(cor(colon_heatmap, brain_heatmap))>.6)
cols_genes_sum_p5 = colSums(abs(cor(colon_heatmap, brain_heatmap))>.5)

gene_colon_p6 =cols_genes_sum[cols_genes_sum==max(cols_genes_sum) ]

cols_genes_sum_p5=colSums(cor(colon_heatmap, brain_heatmap)>.5)
cor_values_p6 =colon_brain_cor_df$ENSG00000138376.6[abs(colon_brain_cor_df$ENSG00000138376.6)>.6]
cor_genes_names_p6 =row.names(colon_brain_cor_df)[abs(colon_brain_cor_df$ENSG00000138376.6)>.6]

ml_genes_df_p6=qn.colon.t[cor_genes_names_p6]
#adding the nrain gene to the colon genes dataframe
ml_genes_df_p6 =cbind(ml_genes_df_p6, qn.brain.t$ENSG00000138376.6[c(0:60)])
##########################