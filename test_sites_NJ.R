library(phangorn)
rm(list=ls())
set.seed(1333343) ## change this for further analyses

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
plot(treeNJ)

init.n <- 150 ## the number of sites

max.tries <- 20000000
max.fails <- max.tries
old.keep <- 0.99
proposals <- matrix(NA, ncol=2, nrow=max.tries)

optimizeSetNJ <- function(data, cog.tree=cog.tree, init.n=ncol(data), max.tries=100000, max.fails=max.tries, minimum.sites=20, maximum.sites=400, old.keep=0.9, proposals=NULL){
    global.datamat=NULL
    global.new.indexes = NULL
    global.curRf = NULL
    ## iterate to get a distribution to have an idea about the
    ## RF distance
    fails <- 0
    #proposals <- matrix(NA, ncol=2, nrow=max.tries)
    ## to get an idea of the random distances
    rand.n = 500
    randomRF = vector("numeric", length=rand.n)
    i <- 1
    for(i in 1:rand.n){
        random.col <- sample(1:ncol(initialdatamat), init.n, replace=TRUE)
        datamat = as.matrix(initialdatamat[,random.col])
        samplephy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)        
        ##optim.pml(treeNJ, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=T, optRate=T)
        dm  <- dist.hamming(samplephy)
        tr <- NJ(dm)
        ##fit = pml(treeNJ, data=samplephy)
        ##tr = optim.pml(fit, model="Mk", rearrangement="stochastic")
        d=treedist(tr, cog.tree, check.labels=TRUE)
        randomRF[i] <- d[1]
    }
    pdf("random_histogram_RF.pdf")
    hist(randomRF) ## gyrw sto 65 mexri peripou to 50
    dev.off()
    cur.indexes = sample(1:ncol(initialdatamat), init.n, replace=FALSE)
    datamat = as.matrix(initialdatamat[,cur.indexes])
    new.indexes = cur.indexes
    i <- 1
    curRf = 1000000
    n.col.to.replace = round(init.n/15, 0)
    i <- 1
    for(i in 1:max.tries){
        if( (i==1) || (i %% 10000 == 0)){
            print(paste("***** d: ", d[1], " Current RF: ", curRf))
            attrNames <- sort(gsub("^X", "", colnames(initialdatamat)[new.indexes]))
            write.table(attrNames, file=paste("selected_features_step", i, ".txt", sep=""), quote=F, col.names=F, row.names=F)
            write.table(curRf, file=paste("current_accuracy_step", i, ".txt", sep=""), quote=F, col.names=F, row.names=F)
            samplephy = as.phyDat(datamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
            dm <- dist.hamming(samplephy)
            tr <- NJ(dm)
            d <- treedist(tr, cog.tree, check.labels =TRUE)
            write.phyDat(x=samplephy, file=paste("current_opt_data_large_step", i, ".nex", sep=""), format="nexus")
            write.table(new.indexes, file=paste("current_opt_indexes_large_step", i, ".csv", sep=""), sep=",", col.names=F, row.names=F, quote=F)
            write.tree(tr, file=paste("current_opt_tree_large_step", i, ".nwk", sep=""))
        }
        ##print(i)
        newdatamat = datamat
        ## how many columns to replace
        ##n.col.to.replace = sample((n.col.to.replace-1):(n.col.to.replace+1), 1)##rbinom(n=1, size=ncol(data), prob=n.col.to.replace/ncol(data))
        ##if(n.col.to.replace == 0){ n.col.to.replace=1 }        
        indexes.to.sample.from <- as.vector(1:ncol(initialdatamat))[-new.indexes]
        ## perc of keeping from the old indexes lets.say oldkeep. This is from a bernoulli process
        old.keep <- max(old.keep, length(new.indexes)/ncol(initialdatamat))
        repeat{
            cur.indexes.to.keep <- sample(c(FALSE, TRUE), size=length(new.indexes), replace=TRUE, prob=c(1-old.keep, old.keep))
            prob1 <- min(length(new.indexes)/length(indexes.to.sample.from)*(1-old.keep), 1)
            new.indexes.to.insert <- sample(c(FALSE,TRUE), size=length(indexes.to.sample.from), prob=c(1-prob1, prob1), replace=TRUE)
            new.indexes.to.propose <- c(new.indexes[cur.indexes.to.keep], indexes.to.sample.from[new.indexes.to.insert])
            if( length(new.indexes.to.propose) > minimum.sites & length(new.indexes.to.propose) < maximum.sites){break}
        }
        newdatamat = as.matrix(initialdatamat)[,new.indexes.to.propose]
        samplephy = as.phyDat(newdatamat, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
        dm <- dist.hamming(samplephy)
        tr <- NJ(dm)
        ##optim.pml(treeNJ, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=T, optRate=T)
        ##fit = pml(treeNJ, data=samplephy)
        ##tr = optim.pml(fit, model="Mk", rearrangement="stochastic", control = pml.control(trace = 0))
        d=treedist(tr, cog.tree, check.labels=TRUE)
        if(d[1] <= curRf || ((d[1] < curRf + 4) && (runif(1, 0, 1) < 0.001))){
            if(d[1] < curRf){
                print(paste("step", i, " found a better matrix... current best: ", curRf, " current size: ", length(new.indexes.to.propose)))
            }
            if(d[1] <= curRf){
                global.datamat=newdatamat
                global.new.indexes = new.indexes.to.propose
                global.curRf = d[1]
                global.tree = tr
                write.tree(tr, file=paste("current_best_tree_large_step", i, ".nwk", sep=""))
            }
            datamat = newdatamat
            new.indexes = new.indexes.to.propose
            curRf = d[1]
            fails = 0
        }else{
            fails  <- fails+1
            #print(c(fails, d[1]))
            #n.col.to.replace = ceiling( round(init.n/10, 0) * (max.fails - fails)/(max.fails+1))
            if(fails > max.fails){
                break
            }
        }
        ##print(paste("***** d: ", d[1], " Current RF: ", curRf, " current size: ", length(new.indexes), " current proposal: ", length(new.indexes.to.propose), "sizes: ", length(new.indexes), ",", sum(cur.indexes.to.keep), ",", sum(new.indexes.to.insert) , sep="") )
        proposals[i,1] <- curRf
        proposals[i,2] <- d[1]
    }
    return(list(datamatrix=datamat, indexes=new.indexes, rf=curRf, size=length(new.indexes), globalmat=global.datamat, global.indexes=global.new.indexes, global.rf = global.curRf))
}

results=optimizeSetNJ(data=initialdatamat, cog.tree=cog.tree, init.n=150, max.tries=max.tries, old.keep=0.99, proposals=proposals)

## get names
attrNames <- gsub("^X", "", colnames(data)[results$indexes])
write.table(attrNames, file="selected_features.txt", quote=F, col.names=F, row.names=F)

samplephy = as.phyDat(results$datamatrix, type='USER', levels=c("0","1"), ambiguity="?", names=rownames(data), n=nrow(data), return.index=TRUE)
dm <- dist.hamming(samplephy)
tr <- NJ(dm)
d <- treedist(tr, cog.tree, check.labels=TRUE)
write.phyDat(x=samplephy, file="opt_data_large.nex", format="nexus")
write.table(results$indexes, file="opt_indexes_large.csv", sep=",", col.names=F, row.names=F, quote=F)
write.tree(tr, file="opt_tree_large.nwk")
d
