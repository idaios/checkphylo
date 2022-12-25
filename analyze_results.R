fileList <- list.files(pattern="current_accuracy_step.*.txt")

minAcc <- 999999999
for(x in fileList){
    a <- read.table(x)
    if(a[1,1]  < minAcc){ minAcc = a[1,1] }
}


highAcc <- sapply(fileList,
                  function(x){
                      a = read.table(x)
                      if(a[1,1] == minAcc){
                          return(1)}
                      return(0)})
    
highAccFiles <- names(highAcc)[highAcc == 1]
inds <- gsub(pattern="[a-z_]+(\\d+).txt", x=highAccFiles, replacement="\\1")
highAccFiles <- names(highAcc)[highAcc == 1][order(as.integer(inds))]
distancesMatrix <- matrix(NA, length(highAccFiles), length(highAccFiles))
colnames(distancesMatrix) <- rownames(distancesMatrix) <- sort(as.integer(inds))

colnames(distancesMatrix)

for(i in 1:nrow(distancesMatrix)){
    for(j in 1:ncol(distancesMatrix)){
        fileA <- colnames(distancesMatrix)[i]
        a <- read.table(paste("selected_features_step", fileA, ".txt", sep=""))[,1]
        fileB <- colnames(distancesMatrix)[j]
        b <- read.table(paste("selected_features_step", fileB,".txt", sep=""))[,1]
        common <- intersect(a,b)
        distancesMatrix[i,j] <- length(a)+length(b) - 2*length(common)
    }
}

library(gplots)

pdf("differences_features_heatmap.pdf")
heatmap.2(distancesMatrix, Rowv=F, Colv=F, trace='none', dendrogram='none')
dev.off()


res <- list()
for(i in 1:length(inds)){
    f <- paste("selected_features_step", inds[i], ".txt", sep="")
    a <- read.table(f)
    res[[i]] <- a[,1]
}

tt <- sort(table(unlist(res)))
tt.top <- tt/length(inds) > 0.0

top.features <- names(tt[tt.top])

top.features

full.data <- read.table("Datasets without geo.xlsx - Full Dataset.csv", sep=",", quote='"')
row.names(full.data) <- full.data[,3]

write.table(cbind(full.data[top.features,c(3,5)], tt[tt.top]/length(inds)), file="best_results00.csv", sep=";", quote=F, row.names=F, col.names=F)

