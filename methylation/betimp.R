###########################################
# betas imputation
###########################################

# B --> beta values
# a450 --> annotations of the array
# clinic --> clinical information of the samples

betimp <- function(B, a450, clinic, resDir){
  
  ###########################################
  # Identify loci with missings and their potential predictors
  nmiss <- apply(B[!rownames(B)%in%a450[a450$chr=='Y', 1], ], 1, FUN=function(o) sum(is.na(o)))
  nmiss.h <- apply(B[rownames(B)%in%a450[a450$chr=='Y', 1], clinic$sex=="Male"], 1, FUN=function(o) sum(is.na(o)))
  psm <- names(nmiss[nmiss > 0])
  psm.h <- names(nmiss.h[nmiss.h >0])
  psmiss <- c(psm, psm.h)
  chrmiss <- a450[psmiss, "chr"]
  posmiss <- a450[psmiss, "pos"]
  
  a450$IlmnID<-rownames(a450)
  
  chru <- unique(chrmiss)
  bps <- c(300, 500, 1000, 2500, 5000)
  ln <-sapply(chru, FUN=function(o)  {
    s <- chrmiss==o
    pm <- posmiss[s]
    names(pm) <- psmiss[s]
    s <- !is.na(a450$chr) & a450$chr==o			
    pos <- a450[s, c("pos")]
    names(pos) <- a450[s, c("IlmnID")]
    pos <- pos[names(pos)%in%rownames(B)]
    dif <- rep(pm, each=length(pos)) - rep(pos, length(pm))
    r <- sapply(bps, FUN=function(o, dif, pos, pm)  {
      wh <- which(abs(dif) < o)
      wha <- wh%%length(pos)
      wha[wha==0] <- length(pos)
      whm <- floor((wh-1)/length(pos)) + 1
      psall <- names(pos[wha])
      psm <- names(pm[whm])
      r <- tapply(psall, psm, FUN=function(o) o)
      r2 <- sapply(1:length(r), FUN=function(j) r[[j]][r[[j]] != names(r)[j]], simplify=F)
      names(r2) <- names(r)
      r2
    }, dif, pos, pm, simplify=F)
    names(r) <- bps
    r
    
  }, simplify=F)
  
  names(ln) <- chru
  
  lneigh <- list()
  for (j in 1:length(ln[[1]]))
  {
    lneigh[[j]] <- list()
    for (k in 1:length(ln))
    {
      lneigh[[j]] <- c(lneigh[[j]], ln[[k]][[j]])
    }
    names(lneigh)[j] <- bps[j]
  }
  
  ###############################################################################################################################
  # Descriptives of number of predictors available at different window sizes
  
  nn <- sapply(lneigh, FUN=function(o) sapply(o, FUN=function(o) length(o)))
  
  sm <- t(apply(nn, 2, FUN=function(o)
  {
    n <- length(o)
    n0 <- length(o[o==0])
    pcn0 <- n0/n*100
    m <- mean(o)
    sd <- sd(o)
    qtc <- quantile(o, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
    c(n=n, n0=n0, pcn0=pcn0, mean=m, sd=sd, qtc)
  }))
  
  sm[, 3:5] <- formatC(sm[, 3:5], format='f', digit=1)
  sm <- cbind(nbp=paste("Nbp =", as.numeric(rownames(sm))*2), sm)
  
  ###############################################################################################################################
  # For probes with less than 5 predictors at 2500bp, select the nearest probes
  
  psadd <- names(which(sapply(lneigh[["2500"]], length) < 5))
  chradd <- a450[psadd, "chr"]
  posadd <- a450[psadd, "pos"]
  
  ln <-sapply(chru, FUN=function(o)
  {
    s <- chradd==o
    pm <- posadd[s]
    names(pm) <- psadd[s]
    s <- !is.na(a450$chr) & a450$chr==o			
    pos <- a450[s, c("pos")]
    names(pos) <- a450[s, c("IlmnID")]
    pos <- pos[names(pos)%in%rownames(B)]
    sapply(pm, FUN=function(o) sort(abs(o - pos))[2:6], simplify=F)
  }, simplify=F)
  names(ln) <- chru
  ladd <- do.call(c, ln)
  names(ladd) <- psadd
  rangeladd <- t(sapply(ladd, range))
  
  smadd <- t(apply(rangeladd, 2, FUN=function(o)
  {
    n <- length(o)
    m <- mean(o)
    sd <- sd(o)
    qtc <- quantile(o, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
    c(n=n, mean=m, sd=sd, qtc)
  }))
  smadd[, -1] <- formatC(smadd[, -1], format='f', digit=1)
  rownames(smadd) <- c("Nearest", "Most distant")
  smadd <- cbind(rownames(smadd), smadd)
  
  ###############################################################################################################################
  # Final list of predictors
  
  ln <- lneigh$"2500"
  ladd.ps <- sapply(ladd, FUN=names, simplify=F)
  ln[names(ladd)] <- ladd.ps
  lneighs.add <- ln
  
  ###############################################################################################################################
  # Predictions with Random Forest
  
  cont <- names(lneighs.add)
  lpreds <- eris.vapply(nomx="cont", nnode=20, FUN=function(o){
    #lpreds <- sapply(cont, FUN=function(o)    {
    ps <- o
    pr <- lneighs.add[[o]]
    chr <- a450[ps, "chr"]
    if (chr == 'Y')
    {
      s <- clinic$sex=='Male'
      
    }else {	
      s <- rep(TRUE, ncol(B))
    }
    
    bs <- B[ps, s]
    
    bpr <- data.frame(t(B[pr, s]), type=clinic$type[s], sex=clinic$sex[s], age=clinic$age[s])
    
    stest <- names(bs[is.na(bs)])
    strain <- names(bs[!is.na(bs)])
    strain <- strain[!strain%in%names(which(apply(bpr, 1, FUN=function(o) any(is.na(o)))))]
    
    bpr.tr <- bpr[strain, ]
    bpr.ts <- bpr[stest, ]
    bs.tr <- bs[strain]
    bs.ts <- bs[stest]
    
    rf <- randomForest(x=bpr.tr, y=bs.tr, ntree=1000, mtry= max(floor(ncol(bpr)/3), 1), importance=T, do.trace=100)
    
    difs <- apply(cbind(bs.tr, rf$predicted), 1, diff)
    mean <- mean(difs)
    sd <- sd(difs)
    qtc <- quantile(difs, c(0, 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99, 1))
    sm <- c(meandif=mean, sddif=sd, meanb=mean(bs.tr), sdb=sd(bs.tr), qtc, mse=rf$mse[length(rf$mse)], rsq=rf$rsq[length(rf$rsq)])
    
    pred <- try(predict(rf, bpr.ts))
    if (class(pred)=='try-error') {
      smiss <- names(apply(bpr.ts, 1, FUN=function(o) any(is.na(o))))
      whnomiss <- sapply(smiss, FUN=function(o) which(!is.na(bpr.ts[o, ])), simplify=F)
      vnomiss <- sapply(whnomiss, FUN=function(o) paste(colnames(bpr.ts)[o], collapse=' + '))
      vnomissu <- unique(vnomiss)
      rfm <- list()
      for (k in 1:length(vnomissu)) {
        vn <- strsplit(vnomissu[k], split=' \\+ ')[[1]]
        rfm[[k]] <- randomForest(x=bpr.tr[, vn], y=bs.tr, ntree=1000, mtry= max(floor(ncol(bpr)/3), 1), importance=T, do.trace=100)
      }
      names(rfm) <- vnomissu
      pred <- rep(NA, length(smiss))
      names(pred) <- smiss
      for (k in 1:length(smiss)) {
        vn <- strsplit(vnomiss[[smiss[k]]], split=' \\+ ')[[1]]
        pred[k] <- predict(rfm[[vnomiss[[smiss[k]]]]], bpr.ts[k, vn])				
      }
      
    }
    
    list(pred=pred, perf=sm, difs=difs)
    
    #}, simplify=FALSE)
  }, simplify=FALSE, obj=c("lneighs.add", "a450", "clinic", "B"), lib=c("randomForest"))
  
  lr <- eris.applies.res(x=lpreds, tjoin='c', remove=F)
  names(lr) <- names(lneighs.add)
  limput <- lr
  
  ###############################################################################################################################
  # Save imputations in the matrix data
  
  preds <- sapply(limput, FUN=function(o) o$pred)
  Betas.imp <- B
  for (j in 1:length(preds))
  {
    print(j)
    Betas.imp[names(preds)[j], names(preds[[j]])] <- preds[[j]]
  }
  
  
  ### Save as Rdata
  save(Betas.imp, file=paste0(resDir, "Betas.imp.Rdata"))
  return(Betas.imp)
  
}
