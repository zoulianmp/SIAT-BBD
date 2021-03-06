# *************课堂作业***********************************
#
# 邹炼 ，医工所微创中心，工号：XJY102908
# lian.zou@siat.ac.cn
#
#**********************************************************



source("./zoulianfuncs/minmaxEuclideanDistance.R")

egDataMatrix <- read.csv(file="class-matrix-GSE10810.csv")

dimension <- dim(egDataMatrix)

dataMatrix <- egDataMatrix[,2:dimension[2]]

samplenames <- dimnames(dataMatrix)[[2]]


  
minmaxindex<-minmaxEuclideanDistance(dataMatrix)

minleft<-minmaxindex[1]
minright<-minmaxindex[2]

maxleft<-minmaxindex[3]
maxright<-minmaxindex[4]

## plot
par(mfrow=c(1,2))
plot( egDataMatrix[,minleft],egDataMatrix[,minright],
     main = "The Minimum of Euclidean Distaces" ,  
     xlab = samplenames[minleft], ylab = samplenames[minright])


abline(0,1)
plot(egDataMatrix[,maxleft], egDataMatrix[,maxright],
     main = "The Maximum of Euclidean Distaces" ,  
     xlab = samplenames[maxleft], ylab = samplenames[maxright])

abline(0,1)

par(mfrow=c(1,1))
