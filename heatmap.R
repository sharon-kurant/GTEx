# brain_heatmap =matrix(as.numeric(t(qn.brain.res[,c(0:60)])),ncol =nrow(qn.brain.res),dimnames=dimnames(t(qn.brain.res[,c(0:60)])))
# colon_heatmap =matrix(as.numeric(t(qn.colon.res)),ncol =nrow(qn.colon.res),dimnames=dimnames(t(qn.colon.res)))
# 
# #colon_brain_heatmap=mapply(cor, as.data.frame(t(qn.colon.res)), as.data.frame(t(qn.brain.res[,c(0:60)])))
# colon_brain_heatmap =cor(colon_heatmap, brain_heatmap)
# heatmap(colon_brain_heatmap)

#tracking which subjects correspond to which columns in the brain dataframe
brainsubtrack=rep(NA, dim(qn.brain.res)[2])
for (i in 1:dim(qn.brain.res)[2]){
  brainsubtrack[i]=strsplit(colnames(qn.brain.res)[i], "-")[[1]][2]
}

#tracking which subjects correspond to which columns in the colon dataframe
colonsubtrack=rep(NaN, dim(qn.colon.res)[2])
for (i in 1:dim(qn.colon.res)[2]){
  colonsubtrack[i]=strsplit(colnames(qn.colon.res)[i], "-")[[1]][2]
}

#tracking which subjects in the brain appear in the colon dataframe
intrack=rep(NA, length(brainsubtrack))
for (i in 1:length(brainsubtrack)){
  intrack[i]=brainsubtrack[i] %in% colonsubtrack
}

#removing subjects which do not appear
qn.brain.res.new=qn.brain.res[,intrack]

#renewing brainsubtracking according to new dataframe with removed columns
brainsubtrack=rep(NA, dim(qn.brain.res.new)[2])
for (i in 1:dim(qn.brain.res.new)[2]){
  brainsubtrack[i]=strsplit(colnames(qn.brain.res.new)[i], "-")[[1]][2]
}

#preparing a colon dataframe to correspond to the brain dataframe
qn.colon.res.braincomp=qn.colon.res[,which(colonsubtrack==brainsubtrack[1])[1]]
for (i in 2:length(brainsubtrack)){
  qn.colon.res.braincomp=cbind(qn.colon.res.braincomp,qn.colon.res[,which(colonsubtrack==brainsubtrack[i])[1]])
}

# n_genes= 5000
# brain_heatmap =matrix(as.numeric(t(qn.brain.res[c(0:5000),c(0:60)])),ncol =5000,dimnames=dimnames(t(qn.brain.res[c(0:5000),c(0:60)])))
# colon_heatmap =matrix(as.numeric(t(qn.colon.res[c(0:5000),])),ncol =5000,dimnames=dimnames(t(qn.colon.res[c(0:5000),])))

#colon_brain_heatmap=mapply(cor, as.data.frame(t(qn.colon.res)), as.data.frame(t(qn.brain.res[,c(0:60)])))
#colon_brain_heatmap =cor(colon_heatmap, brain_heatmap)

#calculating correlation matrix
colon_brain_heatmap=cor(t(qn.colon.res.braincomp), t(qn.brain.res.new))
#heatmap(colon_brain_heatmap)



# colon_brain_cor_df = as.data.frame(colon_brain_heatmap)
# which(colon_brain_cor_df > 0.45)
# 
# colon_brain_cor_df[8,20]
# 
# colnames(colon_brain_heatmap)[apply(colon_brain_heatmap, 1, function (x) which(x==max(x[x>0.5])))]

#calculating the number of correlations whose absolute value is greater than 4 for each gene in brain sample
cols_genes_sum =colSums(abs(colon_brain_heatmap)>.4)
#sorting
cols_genes_sum_sorted=sort(cols_genes_sum, decreasing=TRUE)



#gene.f[which(cols_genes_sum==max(cols_genes_sum)),"Description"]

#genename=rep(NA, length(cols_genes_sum))


#building gene rank dataframe according to the genes which correlate the most to colon samples
generankdf=data.frame(rep(NA, length(cols_genes_sum)))
colnames(generankdf)=c("rank")
generankdf$genedescription=rep(NA, length(cols_genes_sum))
generankdf$numcor=NA

for (i in 1:length(cols_genes_sum)){
  generankdf$rank[i]=i
  generankdf$genedescription[i]=gene.f[,"Description"][gene.f[,"Name"]==row.names(qn.brain.res.new)[order(cols_genes_sum, decreasing=TRUE)[i]]]
  generankdf$genename[i]=row.names(qn.brain.res.new)[order(cols_genes_sum, decreasing=TRUE)[i]]
#  generankdf$genedescription[i]=gene.f[order(cols_genes_sum, decreasing=TRUE)[i],"Description"]
#  generankdf$genename[i]=gene.f[order(cols_genes_sum, decreasing=TRUE)[i],"Name"]
}
generankdf$numcor=cols_genes_sum_sorted


yframe=t(qn.brain.res.new[generankdf$genename[1:10],])

modelX=t(qn.colon.res.braincomp)

connectgeneslists=list()

for (i in 1:10){
  connectgeneslists[[i]]=c(row.names(qn.colon.res.braincomp)[abs(colon_brain_heatmap[,generankdf$genename[i]])>.4])
}
#firstone=row.names(qn.colon.res.braincomp)[abs(colon_brain_heatmap[,generankdf$genename[1]])>.4]

# ylist=data.frame()
# 
# for (i in 1:10){
#   ylist=append(ylist, qn.brain.res.new[generankdf$genename[i],])
# }
# 
# gene_colon_p6 =cols_genes_sum[cols_genes_sum==max(cols_genes_sum) ]
# 
# cols_genes_sum_p5=colSums(cor(colon_heatmap, brain_heatmap)>.5)
# cor_values_p6 =colon_brain_cor_df$ENSG00000138376.6[abs(colon_brain_cor_df$ENSG00000138376.6)>.6]
# cor_genes_names_p6 =row.names(colon_brain_cor_df)[abs(colon_brain_cor_df$ENSG00000138376.6)>.6]
#  
# ml_genes_df_p6=qn.colon.t[cor_genes_names_p6]
# #adding the nrain gene to the colon genes dataframe
# ml_genes_df_p6 =cbind(ml_genes_df_p6, qn.brain.t$ENSG00000138376.6[c(0:60)])
# ##########################
# gene.f[which(gene.f["Name"] == "ENSG00000138376.6")]
# brain_gene = gene.f[gene.f[,"Name"] =="ENSG00000138376.6"]
# colon_genes =gene.f[gene.f[,"Name"] ==cor_genes_names_p6]

