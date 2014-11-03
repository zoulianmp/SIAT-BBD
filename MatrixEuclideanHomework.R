source("./zoulianfuncs/minmaxEuclideanDistance.R")

egDataMatrix <- read.csv(file="class-matrix-GSE10810.csv")

dimension <- dim(egDataMatrix)

dataMatrix <- egDataMatrix[,2:dimension[2]]


minmaxindex<-minmaxEuclideanDistance(dataMatrix)

minleft<-minmaxindex[1]
minright<-minmaxindex[2]

maxleft<-minmaxindex[3]
maxright<-minmaxindex[4]

## plot
par(mfrow=c(1,2))
plot(egDataMatrix[,minleft], egDataMatrix[,minright])
abline(0,1)
plot(egDataMatrix[,maxleft], egDataMatrix[,maxright])
abline(0,1)

par(mfrow=c(1,1))
