library(phangorn)
rm(list=ls())
set.seed(928293) ## change this for further analyses                                                                          

cog.tree=read.tree("IE2011_Cognates_rel_ANNOT.nwk")
data = read.table("46glossesfull.csv", header=TRUE, row.names=1, colClasses='character')
cog.tree$tip.label[cog.tree$tip.label == "Ptg-E"] = "PtgE"
datamat = as.matrix(data)

tab <- apply(datamat, 2, function(x){table(factor(x, levels=c(0,1)))})
indsToUse <- apply(tab, 2, function(x){x[1] > 1 && x[2] > 1})
datamat <- datamat[,indsToUse]
initialdatamat <- datamat

phy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)

dm  <- dist.hamming(phy)
treeNJ <- NJ(dm)

fit <- pml(treeNJ, data=phy)
fitJC <- optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE,rearrangement ="NNI")
plot(density(fitJC$siteLik), col="red", lwd=4)
fit$logLik
fitJC$logLik
siteLikes=c()
plot(density(siteLikes))
for(i in 1:100){
    d = simSeq(fitJC, l=219)
    dm  <- dist.hamming(d)
    tnj <- NJ(dm)
    fit <- pml(tnj, data=d)
    tr = optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE,rearrangement = "NNI")
    siteLikes = c(siteLikes, tr$siteLik)
    points(density(tr$siteLik), col="gray", type='l')
}



plot(fitJC$tree)
names(fit)
