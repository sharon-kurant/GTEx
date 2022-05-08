sampleTree = hclust(dist(t(mat.colon)), method = "average")
# sizeGrWindow(12,9)
par(cex = 0.3)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Colon Clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# calculate cut average for k = 2
cut_avg <- cutree(sampleTree, k = 2)
plot(sampleTree)
rect.hclust(sampleTree , k = 2, border = 2:6)

# Look at cluster results
ind1 = which(cut_avg == 1)
ind2 = which(cut_avg == 2)

colon.clust1.indxs = which(pheno.colon$SAMPID %in% names(ind1))
filtered.pheno.colon = pheno.colon[colon.clust1.indxs,] # 62
filtered.mat.colon = mat.colon[,colon.clust1.indxs]

#filtered.samples.colon = filtered.pheno.colon$SAMPID
#filtered.mat.colon = mat.colon[,colnames(mat.colon) %in% filtered.samples.colon]

colon.clust2.indxs = which(pheno.colon$SAMPID %in% names(ind2))
out.pheno.colon = pheno.colon[colon.clust2.indxs,] # 4

which(pheno.colon$SAMPID %in% out.pheno.colon$SAMPID)
out.subjid.colon = as.vector(out.pheno.colon$SUBJID)

colon.subjid.to.save.idx = which(filtered.pheno.colon$SUBJID %in% out.subjid.colon)
colon.subjid.to.save = filtered.pheno.colon[colon.subjid.to.save.idx,]$SUBJID 

colon.subjid.to.remove.idx = which(!out.subjid.colon %in% colon.subjid.to.save)

remove.this.sample.idx = which(filtered.pheno.brain$SUBJID %in%
                                 out.subjid.colon[colon.subjid.to.remove.idx])
filtered.pheno.brain = filtered.pheno.brain[-remove.this.sample.idx,]
filtered.mat.brain = filtered.mat.brain[,-remove.this.sample.idx]

colon.to.remove.idx = which(!unique(filtered.pheno.colon$SUBJID) %in% unique(filtered.pheno.brain$SUBJID))
filtered.pheno.colon = filtered.pheno.colon[-colon.to.remove.idx,]
filtered.mat.colon = filtered.mat.colon[,-colon.to.remove.idx]
