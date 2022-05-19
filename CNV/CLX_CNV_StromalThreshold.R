
## -------------------------------------------------------------------------------------------

## library
library(gdata)
library(beeswarm)

## -------------------------------------------------------------------------------------------

## Stroma
dat.samp = read.xls("/shared/projects/COLONOMICS/data/CLX_Samples.xls",head=T)
dat.samp = subset(dat.samp, type=="Tumor")
rownames(dat.samp) = dat.samp$id_clx_individual

## -------------------------------------------------------------------------------------------

strome = dat.samp$stromal_score
names(strome) = rownames(dat.samp)
strome = strome[!(is.na(strome))]

## Histograma
hist(strome, prob=T)
lines(density(strome))

## Cluester
d = dist(strome, method = "euclidean") 
fit = hclust(d, method="ward.D") 

## 2 grupos
par(mfrow=c(1,2))
plot(fit)
groups = cutree(fit, k=2)
rect.hclust(fit, k=2, border="red")
boxplot(strome ~  groups)
abline(h=0)

## 3 grupos
par(mfrow=c(1,2))
plot(fit)
groups = cutree(fit, k=3)
rect.hclust(fit, k=3, border="red")
boxplot(strome ~  groups)
abline(h=0)
abline(h=-1250)

## 4 grupos
groups = cutree(fit, k=4)
tiff("hc.tiff", width = 2000,height = 2000,unit='px', res=300)	
par(mar=c(2,4,2,2))
plot(fit, main="Cluster dendrogram CLX", cex=0.5, xlab="",sub="")
rect.hclust(fit, k=4, border="red")
dev.off()
groups.cof = ifelse(groups==1, 0.2, ifelse(groups==2,0.5, ifelse(groups==3,0.4,0.3)))
tiff(file="Supplementary Figure 1.tiff", width = 2000,height = 2000,unit='px', res=300)	
par(oma=c(0,0,0,0), mar=c(5,5,2,2))
boxplot(strome ~  groups.cof, xlab="Cutoff of CNE", ylab="Proportion Stroma", lwd=2, cex.lab=1.5)
title("Cutt-off value of CLX")
abline(h=0, lty=2, col="grey15", lwd=2)
abline(h=-1250, lty=2, col="grey15", lwd=2)
abline(h=800, lty=2, col="grey15", lwd=2)
dev.off()

## -------------------------------------------------------------------------------------------

dat.samp$strome = ifelse(dat.samp$stromal_score<(-1250),0.5,ifelse(dat.samp$stromal_score<0,0.4,ifelse(dat.samp$stromal_score<800, 0.3,0.2)))
dat.samp$strome[is.na(dat.samp$strome)] = 0.4
p.cut = data.frame(ID=dat.samp$id_clx_individual, p=dat.samp$strome)
rownames(p.cut) = p.cut$ID
save(p.cut, file="p.cut.RData")