#############################################################
# COLONOMICS METHYLATION DATA PREPROCESSING
#############################################################
# JUNE 2013
#############################################################


library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(gdata)


#############################
# Directories

resDir <- "/mnt/hydra/ubs/shared/users/Anna/COLONOMICS/MultiOmics/PreprocessCode/"
baseDir <- "/mnt/hydra/ubs/shared/projects/COLONOMICS/data/methylation/ICO/Intensity_Data_Files/"
PriceDir <- "/mnt/hydra/ubs/shared/projects/COLONOMICS/methylation/data/Annotation/"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Normalize Betas:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################
# Remove samples not passing QC

# samples data
targets <- read.metharray.sheet(baseDir, recursive=TRUE, verbose=TRUE)

# samples to remove
samp2rm <- c("E2138_M","F2123_M")
targets <- targets[!targets$Sample_Name%in%samp2rm,]

# remove bloods
targets <- targets[-grep("_B",targets$Sample_Name),]

# fix ids
targets$Sample_Name[targets$Sample_Name=="Z2015_N_0"] <- "Z2015_N"
targets$Sample_Name[targets$Sample_Name=="Z2015_T_0"] <- "Z2015_T"

#############################
# Read raw data
RGset <- read.metharray.exp(base=NULL, targets=targets, extended=TRUE, recursive=TRUE)

# samples info
pd <- pData(RGset)
pd$type <- factor(substr(pd$Sample_Name,7,7), levels=c("T","M","N"))
pd$ind <- substr(as.character(pd$Sample_Name),1,5)

# QC report
# qcReport(RGset, sampNames=pd$Sample_Name, sampGroups=pd$Sample_Plate, pdf=paste0(resDir,"qcReport.pdf"))

#############################
# Subset-quantile Within Array Normalization (SWAN)

MSet <- preprocessRaw(RGset)
MSet <- preprocessSWAN(RGset, mSet=MSet)

# Get betas
Betas<-getBeta(MSet)

#############################
# Set as missing if detection pval > 0.01
pval<-detectionP(RGset)>0.01
Betas[pval]<-NA

# Filter out probes with more than 5% of detection pval > 0.01
nmiss<-apply(pval,1,sum)
psrm<-names(nmiss[nmiss>ncol(pval)*0.05])
Betas<-Betas[!rownames(Betas)%in%psrm, ]
rm(MSet,nmiss,pd,psrm,pval)
gc()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Annotations:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#############################
# Remove SNPs
RGset2 <- mapToGenome(RGset)
RGset2 <- addSnpInfo(RGset2)
RGset2 <- dropLociWithSnps(RGset2, snps=c("SBE","CpG"), maf=0)

Betas<-Betas[rownames(Betas)%in%rownames(RGset2),]

rm(RGset2)
gc()

#############################
# remove CpGs mapping to multiple locations

# Multiple alignment of probes -> from E Magda Price et al, Epigenetics and Chromatine 2013: Additional annotation 
# enhances potential for biologically-relevant analysis of the Illumina Infinium HumanMethylation450 BeadChip array
# Annotation downloaded from GEO in 2013/06/13 -> platform GPL16304, 
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304 

# Read enhanced annotation file
price <- read.table(paste0(PriceDir,"GPL16304-47833.txt"), sep='\t', header=T, as.is=T, comment.char='#')

#############################
# Check and save number of probes showing multi-aligment

price <- price[!((!is.na(price$AlleleA_Hits) & price$AlleleA_Hits > 1) | (!is.na(price$AlleleB_Hits) & price$AlleleB_Hits > 1)),]

#############################
# Remove CpHs and controls

price <- price[substring(price$ID,1,2)%in%c("cg"),]

##############################
# merge with locations

data(Locations)
annot <- merge(Locations, price, by.x="row.names", by.y="ID", all=FALSE)
# 441142

rownames(annot) <- annot$Row.names
annot$Row.names <- NULL

rm(price,Locations)
gc()


# Filter betas not in annot
Betas <- Betas[rownames(Betas)%in%rownames(annot),]
annot <- annot[rownames(annot)%in%rownames(Betas),]
Betas <- Betas[rownames(annot),]

colnames(Betas) <- gsub("CLX_Methylation_","",colnames(Betas))
colnames(Betas)[colnames(Betas)=="Z2015_N_0"] <- "Z2015_N"
colnames(Betas)[colnames(Betas)=="Z2015_T_0"] <- "Z2015_T"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Impute missings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

clinic <- read.xls("/mnt/hydra/ubs/shared/projects/COLONOMICS/data/CLX_Samples.xls")
rownames(clinic) <- clinic[,1]
clinic <- clinic[colnames(Betas),]
clinic$type <- factor(substr(rownames(clinic),7,7), levels=c("T","M","N"))

source("betimp.R")
betimp(B=Betas, a450=annot, clinic=clinic, resDir=resDir)
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sexual chromosomes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Chomosome Y in females to missing
Betas.imp[rownames(annot)[annot$chr=='Y'], rownames(clinic)[clinic$sex=='Female']] <- NA

save(Betas.imp, file=paste0(resDir,"betas.Rdata"))



