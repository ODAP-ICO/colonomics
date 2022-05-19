############################################################
# COLONOMICS COPY NUMBER VARIATION
############################################################
# JUNE 2015
############################################################

library("Vega")
library(parallel)
cores <- 30 
library(data.table)
library(gdata)

dir_dat="/shared/projects/COLONOMICS/cn/userHA/APT/"

## -------------------------------------------------------------------------------------------

## merge consecutive regions without CNV 
clear_seg = function(segm, lrr, p.cut){
  
  lrr$Position = as.numeric(paste0(lrr$Position))
  lrr[,4] = as.numeric(paste0(lrr[,4]))
  
  reg = 1
  changes = 1
  for(i in 2:NROW(segm)){
    if(segm[i,"Label"]!=0 | segm[i,"Label"]!=segm[i-1,"Label"]){ reg = reg+1}
    changes= c(changes, reg) 
  }
  
  ## new seg
  seg.new = data.frame(Start=aggregate(segm[,2], by=list(changes), min)$x, End=aggregate(segm[,3], by=list(changes), max)$x)
  seg.new = data.frame(Chromosome=rep(segm[1,1],NROW(seg.new)), seg.new)
  seg.new$Num_Markers = sapply(1:NROW(seg.new), function(x) sum(lrr$Position>=seg.new[x,2] & lrr$Position<=seg.new[x,3]))
  seg.new$Mean = sapply(1:NROW(seg.new), function(x) mean(lrr[lrr$Position>=seg.new[x,2] & lrr$Position<=seg.new[x,3],4]))
  seg.new$Label = ifelse(seg.new$Mean<(-p.cut), -1, ifelse(seg.new$Mean>p.cut, 1, 0))
  
  return(seg.new)
}


## region evaluation
join_seg = function(segm, lrr, p.cut){
  
  join = sapply(2:(NROW(segm)), function(i){
    lrr.1 = lrr[lrr$Position>=segm$Start[i-1] & lrr$Position<=segm$End[i-1],4]
    lrr.2 = lrr[lrr$Position>=segm[i,2] & lrr$Position<=segm[i,3],4]
    if(abs(mean(lrr.1)-mean(lrr.2))>0.5){FALSE}else{
      ifelse(t.test(lrr.1, lrr.2)$p.value<0.0001, FALSE,TRUE)}
  })
  
  #join = c(TRUE, join)
  reg = 1
  changes = 1
  for(i in join){
    if(i==F){ reg = reg+1}
    changes= c(changes, reg) 
  }
  
  seg.new = data.frame(Start=aggregate(segm$Start, by=list(changes), min)$x, End=aggregate(segm$End, by=list(changes), max)$x)
  seg.new = data.frame(Chromosome=rep(segm[1,1],NROW(seg.new)), seg.new)
  seg.new$Num_Markers = sapply(1:NROW(seg.new), function(x) sum(lrr$Position>=seg.new$Start[x] & lrr$Position<=seg.new$End[x]))
  seg.new$Mean = sapply(1:NROW(seg.new), function(x) mean(lrr[lrr$Position>=seg.new$Start[x] & lrr$Position<=seg.new$End[x],4]))
  seg.new$Label = ifelse(seg.new$Mean<(-p.cut), -1, ifelse(seg.new$Mean>p.cut, 1, 0))
  
  return(seg.new)
}

## -------------------------------------------------------------------------------------------

## data cnv
### Tumor
files <- Sys.glob("res_ind/T/*")
listOfFiles <- lapply(files, function(x) fread(x))

tumor.LRR <- do.call(cbind, sapply(1:NROW(listOfFiles), function(x) as.vector(listOfFiles[[x]][,5, with=F])))
tumor.LRR <- data.frame(data.frame(listOfFiles[[1]][,1:3, with=F]), tumor.LRR)
colnames(tumor.LRR) <- c("Name","Chromosome","Position",substring(files, 11, 15))
ind <- colnames(tumor.LRR)[-c(1:3)]

## centromere
centromere = c(124300000,93300000,91700000,50700000,47700000,60500000,59100000,45200000,51800000,40300000,52900000,
               35400000,0,0,0,38200000,22200000,16100000,28500000,27100000,0,0,59500000,13000000)

## Stromal threshold
load("p.cut.RData")

## -------------------------------------------------------------------------------------------

## Normalization (Smooth)
tumor.LRR.2 = tumor.LRR[,c(1:3)]
for(ii in 4:NCOL(tumor.LRR)){
  dat.ind = tumor.LRR[,c(1:3,ii)]
  res = lapply(1:22, function(chr){
    dat.chr = subset(dat.ind, Chromosome==chr)
    dat.smooth = smooth(dat.chr[,4])
    as.numeric(dat.smooth)
  })
  res = unlist(res)
  tumor.LRR.2 = data.frame(tumor.LRR.2, res)
}
colnames(tumor.LRR.2) = colnames(tumor.LRR)
tumor.LRR = tumor.LRR.2
save(tumor.LRR, file="tumor.LRR_smooth.RData")
rm(tumor.LRR.2)

## Stoma
table(colnames(tumor.LRR)[-(1:3)] %in% p.cut$ID)
p.cut = p.cut[colnames(tumor.LRR)[-(1:3)], ]
table(colnames(tumor.LRR)[-(1:3)] == p.cut$ID)

## segmentation
fun = function(j){
  dat.aux = tumor.LRR[,c(2,3,3,j)]
  all.seg = NULL
  for(chr in 1:22){
    
    ## Extract chr data
    dat.aux.chr = dat.aux[dat.aux[,1]==as.character(chr),] 
    dat.aux.chr = dat.aux.chr[order(as.numeric(dat.aux.chr[,2])),]
    
    ## Centromere
    dat.p = dat.aux.chr[as.numeric(dat.aux.chr[,2])<=centromere[chr],]
    dat.q = dat.aux.chr[as.numeric(dat.aux.chr[,2])>centromere[chr],]
    
    ## Segmentation
    seg = vega(CNVdata=dat.q, chromosomes=c(chr), thr_seg=c(p.cut$p[j-3]), beta=0.7)
    colnames(seg) = c("Chromosome","Start","End","Num_Markers","Mean","Label")
    seg.aux = 1
    while(NROW(seg)!=NROW(seg.aux) & NROW(seg)>1){
      seg = clear_seg(seg, dat.q, p.cut$p[j-3])
      seg.aux = seg
      if(NROW(seg)>1){seg = join_seg(seg.aux, dat.q, p.cut$p[j-3])}
    }
    if(NROW(dat.p)!=0){
      seg.p = vega(CNVdata=dat.p, chromosomes=c(chr), thr_seg=p.cut$p[j-3], beta=0.7)
      colnames(seg.p) = c("Chromosome","Start","End","Num_Markers","Mean","Label")
      seg.aux = 1
      while(NROW(seg.p)!=NROW(seg.aux) & NROW(seg.p)>1){
        seg.p = clear_seg(seg.p, dat.p, p.cut$p[j-3])
        seg.aux = seg.p
        if(NROW(seg.p)>1){seg.p = join_seg(seg.aux, dat.p, p.cut$p[j-3])}
      }
      seg = rbind(seg.p, seg)
    }
    
    ## results chr		
    all.seg = rbind(all.seg, seg)
    print(paste0("chr ",chr))
  }
  cbind(rep(ind[j-3], NROW(all.seg)), all.seg)
}

segment = mclapply(X=4:NCOL(tumor.LRR), FUN=fun, mc.cores=cores)
segment = rbindlist(segment)
colnames(segment)[1] = "ID"
for(i in 2:6) segment[,i] = as.numeric(paste0(segment[,i]))

save(segment, file="CLX_CNV_segment.RData")


