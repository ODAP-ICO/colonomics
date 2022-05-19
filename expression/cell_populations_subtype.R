###-----------------------------------------------------------------------------
### Rebeca S.P.
### February 2017
### ESTIMATE stromal content
###-----------------------------------------------------------------------------

library(estimate)

# load CLX data
load("CLX_Expression_PC1.RData")
write.table(exp.genes, "estimate.txt", sep="\t", col.names=TRUE, row.names=TRUE,quote=FALSE)

filterCommonGenes("estimate.txt", output.f="estimate_output.gct", id="GeneSymbol")
estimateScore("estimate_output.gct", "endo_estimate_score.gct", platform="affymetrix")


###-----------------------------------------------------------------------------
### Susanna Aussó
### February 2017
### ESTIMATE stromal content
###-----------------------------------------------------------------------------

