# brainTissues = c("Brain - Amygdala",	"Brain - Anterior cingulate cortex (BA24)", 
#                  "Brain - Caudate (basal ganglia)",	"Brain - Cerebellar Hemisphere",	
#                  "Brain - Cerebellum",	"Brain - Cortex",	"Brain - Frontal Cortex (BA9)",
#                  "Brain - Hippocampus", "Brain - Hypothalamus",	"Brain - Nucleus accumbens (basal ganglia)",
#                  "Brain - Putamen (basal ganglia)",	"Brain - Substantia nigra",	
#                  "Brain - Spinal cord (cervical c-1)",	"Brain - Substantia nigra")

brainTissues = c("Brain - Cerebellum",	"Brain - Cortex",
                 "Brain - Hippocampus", "Brain - Hypothalamus")

bTissue = c("Brain - Cerebellum")

brainData <- pheno.f[pheno.f$'SMTSD' %in% brainTissues,] #1256 rows
colonData <- pheno.f[pheno.f['SMTS'] == "Colon",] #196 rows

x <- c()
y <- c()

subjID <- "GTEX-12ZZX"

sampBrain <- as.vector(brainData[brainData['SUBJID'] == subjID,]$SAMPID) # 96 sampleID's
sampColon <- as.vector(colonData[colonData['SUBJID'] == subjID,]$SAMPID) # 196 sampleID's


par(mfrow=c(2,2))
for (i in 1:4){
  bTissue = pheno.f[pheno.f$'SAMPID' == sampBrain[i], 'SMTSD']
  plot(mat.f.coding[,sampColon], mat.f.coding[,sampBrain[i]], xlab =  "Colon",
       ylab =  bTissue, main = paste(bTissue," vs. ","Colon"))
}

pheno.f[pheno.f$'SAMPID' == sampBrain[1], 'SMTSD']

