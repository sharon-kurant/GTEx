load("gene.f.RData")
load("gene.f.with.entrez.RData")
load("mat.f.coding.RData")
load("pheno.f.RData")

pheno.brain <- pheno.f[pheno.f['SMTS'] == "Brain",] # 1259
pheno.colon <- pheno.f[pheno.f['SMTS'] == "Colon",] # 345

samplesBrain <- c()
samplesColon <- c()
SubjId <- c()

sampBrain <- as.vector(pheno.brain$SAMPID) # 1259 sampleID's
sampColon <- as.vector(pheno.colon$SAMPID) # 345 sampleID's

subjBrain <- as.vector(pheno.brain$SUBJID) # 1259 sujectID's
subjColon <- as.vector(pheno.colon$SUBJID) # 345 subjectID's

for (i in 1:nrow(pheno.brain)) {
  for (j in 1:nrow(pheno.colon)) {
    if (subjBrain[i] == subjColon[j]) {
      samplesBrain <- append(samplesBrain, sampBrain[i])
      samplesColon <- append(samplesColon, sampColon[j])
      SubjId <- append(SubjId, subjColon[j])
    }
  }
}

SubjUnique <- unique(SubjId)
bcSampsID <- append(unique(samplesBrain), unique(samplesColon))
  
pheno.bc <- pheno.f[pheno.f$SAMPID %in% bcSampsID,]
pheno.brain <- pheno.brain[pheno.brain$SAMPID %in% unique(samplesBrain),] # 353
pheno.colon <- pheno.colon[pheno.colon$SAMPID %in% unique(samplesColon),] # 66

gc()
memory.limit(10000)

mat.bc = mat.f.coding[,bcSampsID]
mat.brain = mat.f.coding[,unique(samplesBrain)]
mat.colon = mat.f.coding[,unique(samplesColon)]

rm(subjBrain, subjColon, mat.f.coding)
gc()

# par(mfrow=c(2,2))
# num = 70
# plot(mat.f.coding[,x[num]], mat.f.coding[,y[num]], xlab =  "colon genes", ylab =  "brain genes", main = subjBrain[num])
# num = 6
# plot(mat.f.coding[,x[num]], mat.f.coding[,y[num]], xlab =  "colon genes", ylab =  "brain genes", main = subjBrain[num])
# num = 150
# plot(mat.f.coding[,x[num]], mat.f.coding[,y[num]], xlab =  "colon genes", ylab =  "brain genes", main = subjBrain[num])
# num = 99
# plot(mat.f.coding[,x[num]], mat.f.coding[,y[num]], xlab =  "colon genes", ylab =  "brain genes", main = subjBrain[num])

