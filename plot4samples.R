brainData <- pheno.f[pheno.f['SMTS'] == "Brain",] #96 rows
colonData <- pheno.f[pheno.f['SMTS'] == "Colon",] #196 rows

x <- c()
y <- c()


sampBrain <- as.vector(brainData$SAMPID) # 96 sampleID's
sampColon <- as.vector(colonData$SAMPID) # 196 sampleID's

subjBrain <- as.vector(brainData$SUBJID) # 96 sujectID's
subjColon <- as.vector(colonData$SUBJID) # 196 subjectID's

for (i in 1:nrow(brainData)) {
  for (j in 1:nrow(colonData)) {
    if (subjBrain[i] == subjColon[j]) {
      x <- append(x, sampBrain[i])
      y <- append(y, sampColon[j])
    }
  }
}


par(mfrow=c(2,2))
num = 70
plot(mat.f.coding[,x[num]], mat.f.coding[,y[num]], xlab =  "colon genes", ylab =  "brain genes", main = subjBrain[num])
num = 6
plot(mat.f.coding[,x[num]], mat.f.coding[,y[num]], xlab =  "colon genes", ylab =  "brain genes", main = subjBrain[num])
num = 150
plot(mat.f.coding[,x[num]], mat.f.coding[,y[num]], xlab =  "colon genes", ylab =  "brain genes", main = subjBrain[num])
num = 99
plot(mat.f.coding[,x[num]], mat.f.coding[,y[num]], xlab =  "colon genes", ylab =  "brain genes", main = subjBrain[num])

