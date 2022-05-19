###############################################################################
#                                                                             #
#  COLONOMICS  -  APRIL 2014                                                  #
#                                                                             #
#  Code for reading and normalizing the COLONOMICS expression arrays.         #
#  And to filter out bad quality samples.                                     #   
#                                                                             #
###############################################################################


# Loading libraries -----------------------------------------------------------
library(affy)
library(limma)
library(mclust)
# -----------------------------------------------------------------------------


# Reading & preparing the clinical data  --------------------------------------
infoClin <- read.table("/shared/projects/COLONOMICS/data/expression/ICO/CLX_RNA_infoClinic.txt", sep="\t", header=TRUE, as.is=TRUE, strip.white=TRUE)
cel.files.list <- c(infoClin$id_hyb_mucosa[infoClin$type=="Mucosa"], infoClin$id_hyb_normal[infoClin$type=="Case"],infoClin$id_hyb_tumor[infoClin$type=="Case"])
cel.files.list <- c(cel.files.list,"CLX_RNA_V2087_N_0.CEL")
cel.files.list <- c(cel.files.list,"CLX_RNA_V2087_N_2.CEL")
cel.files.list <- c(cel.files.list,"CLX_RNA_V2087_T_0.CEL")
cel.files.list <- c(cel.files.list,"CLX_RNA_V2087_T_2.CEL")
infoClin$id_hyb_mucosa <- sub("CLX_RNA_","",sub(".CEL","",infoClin$id_hyb_mucosa))
infoClin$id_hyb_normal <- sub("_1","",sub("CLX_RNA_","",sub(".CEL","",infoClin$id_hyb_normal)))
infoClin$id_hyb_tumor <- sub("_1","",sub("CLX_RNA_","",sub(".CEL","",infoClin$id_hyb_tumor)))

# -----------------------------------------------------------------------------

# Reading sample data  ----------------------------------------
infoSamples <- read.table("/shared/projects/COLONOMICS/data/expression/ICO/CLX_RNA_infoSamples.txt", sep="\t", header=TRUE, as.is=TRUE, strip.white=TRUE)

# -----------------------------------------------------------------------------

# RMA Normailzation ---------------------
raw.data <- ReadAffy( filenames=paste("/shared/projects/COLONOMICS/data/expression/ICO/CEL_files/", cel.files.list, sep="") )
rma.data <- rma(raw.data)
exp.data <- exprs(rma.data)
exp.data <- round(exp.data,6)


# -----------------------------------------------------------------------------

# Assign short colnames -----------------------
exp.data <- exp.data[ , !colnames(exp.data) %in% cel.files.list[grep("_(0|2).CEL",cel.files.list)] ]
colnames(exp.data) <- sub("_1","",sub("CLX_RNA_","",sub(".CEL","",colnames(exp.data))))

# -----------------------------------------------------------------------------

# Reading U219 Array annotation (Release 34 - 10/24/13)) ---------------------
annot <- read.csv("HG-U219.na34.annot.csv", header=TRUE, as.is=TRUE, strip.white=TRUE, skip=25)
rownames(annot) <- annot$Probe.Set.ID

# Set the annotation probesets with the same order that expression data.
annot <- annot[rownames(exp.data),]

genes<-strsplit(annot$Gene.Symbol, " /// ")
genes.unique<- unique(unlist(genes))
names(genes)<-rownames(annot)
gene.probes<-lapply(1:length(genes), function(i) {
  x<- data.frame(genes=genes[[i]])
  x$probe <- names(genes[i])
  x
})
genes.probes<-do.call("rbind",gene.probes)


# -----------------------------------------------------------------------------

# Remove bad quality samples ---------------
samplesForRemove <- rownames(infoSamples[infoSamples$qc_summary=="BAD",])
exp.data <- exp.data[, !colnames(exp.data) %in% samplesForRemove ]

# -----------------------------------------------------------------------------

# Save probeset expression data
save("exp.data",file=paste0("CLX_Expression.RData"))

# -------------------------------------------------------------------------------------------------------------------

# Summary at gene level


# Reading U219 Affy Array annotation (Release 34 - 10/24/13)) ---------------------
annot <- read.csv(gzfile("../annot/HG-U219.na34.annot.csv.gz"), header=TRUE, as.is=TRUE, strip.white=TRUE, skip=25)
rownames(annot) <- annot$Probe.Set.ID

# Set the annotation probesets with the same order that expression data.
annot <- annot[rownames(exp.data),]


library(hgu219.db)

mapped_probes <- mappedkeys(hgu219SYMBOL) # Get the probe identifiers that are mapped to a gene symbol

gene.symbol <- unlist(lapply(as.list(hgu219SYMBOL[mapped_probes]),function(x)x[[1]]))
gene.chr <- as.character(factor(unlist(lapply(as.list(hgu219CHR[mapped_probes]),function(x)x[[1]]))))
gene.ini <- unlist(lapply(as.list(hgu219CHRLOC[mapped_probes]),function(x)x[[1]]))
gene.end <- unlist(lapply(as.list(hgu219CHRLOCEND[mapped_probes]),function(x)x[[1]]))
gene.map <- unlist(lapply(as.list(hgu219MAP[mapped_probes]),function(x)x[[1]]))
gene.genename <- unlist(lapply(as.list(hgu219GENENAME[mapped_probes]),function(x)x[[1]]))
gene.ensembl <- unlist(lapply(as.list(hgu219ENSEMBL[mapped_probes]),function(x)x[[1]]))
gene.entrez <- unlist(lapply(as.list(hgu219ENTREZID[mapped_probes]),function(x)x[[1]]))
gene.strand <- sign(gene.ini)

annot.db <- data.frame(gene.symbol, probe=mapped_probes, chr=gene.chr, start=abs(gene.ini), end=abs(gene.end), strand=gene.strand, map=gene.map, name=gene.genename,  ensembl=gene.ensembl, entrez=gene.entrez  )

genes.symbol<-annot$Gene.Symbol
names(genes.symbol)<-rownames(annot)
gene.symbol.unique<-unique(gene.symbol)

# mapped
genes.symbol.mapped<-genes.symbol[mapped_probes]

# nomapped
genes.symbol.nomapped<-genes.symbol[nomapped]


# recover by gene symbol

library(biomaRt)

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

genes.nomapped <- strsplit(genes.symbol[nomapped], " /// ")
nomapped.genes <- lapply(1:length(genes.nomapped), function(i) {
  x <- data.frame(gene.symbol=genes.nomapped[[i]])
  x$probe <- names(genes.nomapped[i])
  x
})
nomapped.genes <- do.call("rbind",nomapped.genes)

nomapped.genes.repe <- nomapped.genes[  nomapped.genes$gene.symbol %in% gene.symbol, ]
nomapped.genes <- nomapped.genes[  nomapped.genes$gene.symbol %nin% gene.symbol, ]


gene.nomap <- getBM(filters="hgnc_symbol", values=unique(nomapped.genes$gene.symbol), attributes=c("hgnc_symbol","chromosome_name", "start_position","end_position", "strand", "band", "description", "ensembl_gene_id",  "entrezgene"), mart=mart)
colnames(gene.nomap)[1] <- "gene.symbol"
gene.nomap <- gene.nomap[gene.nomap$chromosome_name %in% c(1:22,"X","Y", "MT"),]
gene.nomap <- gene.nomap[!duplicated(gene.nomap$gene.symbol),]

probes.nomapped <- merge( nomapped.genes, gene.nomap, by="gene.symbol") # multiple genes/probe

# fix band & description
probes.nomapped$band <- paste0(probes.nomapped$chromosome_name, probes.nomapped$band )

desc <- strsplit(probes.nomapped$description," [", fixed=TRUE) 
desc <- unlist( lapply(desc, function(x)x[[1]]))
probes.nomapped$description <- desc

probes.orphan <- rownames(annot)[ rownames(annot) %nin% c( mapped_probes,unique(probes.nomapped$probe), unique(nomapped.genes.repe$probe) ) ]
names(probes.nomapped) <- names(annot.db)

# combine
genes.annot <- rbind(annot.db, probes.nomapped)
genes.annot <- genes.annot[!duplicated(genes.annot$gene.symbol),]
rownames(genes.annot) <- genes.annot$gene.symbol
genes.annot$probe <- NULL

# PC1 expression
genes.probes <- rbind(annot.db[,1:2], probes.nomapped[,1:2])
exp.genes <- lapply( genes.annot$gene.symbol, function(x){
  m <- exp.data[ genes.probes[genes.probes$gene.symbol==x,2] , ,drop=FALSE]
  if(nrow(m)>1){
    maxsd <- max(apply(m,1,sd))
    p <- prcomp(t(m)) 
    w <- p$rotation[,1]^2 
    x <- p$x[,1]
    x <- x*sign(cor(x, apply(m,2,mean)))
    m <- x+sum(apply(m,1,mean)*w)
    s <- min(sd(m),maxsd)
    m <- (m-mean(m))*s/sd(m) + mean(m) # adjust sd to sd of probe with max sd if < sd(pc1)
    if(min(m)<0){ # if min is negative, adjust more sd
      d <- min(m)
      m <- (m-mean(m))*mean(m)/(mean(m)-d) + mean(m)
      cat(" - ", x, "\n")
    }
  }
  m
})

exp.genes <- do.call(rbind,exp.genes)
rownames(exp.genes) <- genes.annot$gene.symbol

save(exp.genes, file="CLX_Expression_PC1.RData")
