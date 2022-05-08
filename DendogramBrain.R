sampleTree = hclust(dist(t(mat.brain)), method = "average")
# sizeGrWindow(12,9)
par(cex = 0.3)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Brain Clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# calculate cut average for k = 2
cut_avg <- cutree(sampleTree, k = 2)
plot(sampleTree)
rect.hclust(sampleTree , k = 2, border = 2:6)

# Look at cluster results
ind1 = which(cut_avg == 1)
ind2 = which(cut_avg == 2)

brain.clust1.indxs = which(pheno.brain$SAMPID %in% names(ind1))
filtered.pheno.brain = pheno.brain[brain.clust1.indxs,] # 285
filtered.mat.brain = mat.brain[,brain.clust1.indxs]

#filtered.samples.brain = filtered.brain.pheno$SAMPID
#filtered.mat.brain = mat.brain[,colnames(mat.brain) %in% filtered.samples.brain]

brain.clust2.indxs = which(pheno.brain$SAMPID %in% names(ind2))
out.pheno.brain = pheno.brain[brain.clust2.indxs,] # 68

which(pheno.brain$SAMPID %in% out.pheno.brain$SAMPID)
