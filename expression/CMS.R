
library(CMScaller)

# Load CLX expression data
load("data/CLX_Expression_PC1.RData")
# 20070   246

# Only Tumor samples:
CLX.T<-exp.genes[, grep("_T", colnames(exp.genes))]
# 20070    98

# CMSclassifier requires Entrez ID as rownames: 
rownames(CLX.T)<-genes.annot[rownames(CLX.T), "entrez"]

CLX.T<-CLX.T[!is.na(rownames(CLX.T)),]

res <- CMScaller(CLX.T, RNAseq=FALSE, doPlot=FALSE)
resCMS<-cbind(resCMS,res[,1])
colnames(resCMS)[5]<-"CMScaller"
save(resCMS,file="CLX_CMS.RData")

