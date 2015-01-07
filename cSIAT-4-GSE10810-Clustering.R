#-------------------------------------------------------------------------------
# R scripts for demonstration of the algorithms taught in the course SIAT-BMI
# Teacher: Fengfeng Zhou
# Email: FengfengZhou@gmail.com
# Web: http://healthinformaticslab.org/course/SIAT-BMI/
# Semester: August 2014
# Version 1.0.1
# Update: 2014-08-25
# Dataset: GSE32175
#   Development of Transcriptomic Biomarker Signature in Human Saliva to Detect
#   Lung Cancer
#   URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32175
# Restriction and disclaimer: 
#   These scripts are provided for the course project only. Please do not
#   re-distribute all these scripts, in part or completeness. For any other
#   purpose (including both commercial and non-profit purposes), please
#   contact the teacher for authorization and license.
#   
#   This section must be retained with all these course scripts.
#   
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#-------------------------------------------------------------------------------

#----------------------------------
#Zou Lian  Modified Codes 
#----------------------------------


# Initialization: loading the dataset
egMatrix <- read.csv("class-matrix-GSE10810.csv", header=TRUE, sep=",", row.names=1);
egClass <- read.csv("class-GSE10810-binary.csv", header=TRUE, sep=",", row.names=1);
indexP <- which(egClass$Class == "P");
indexN <- which(egClass$Class == "N");

#egNumSmall <- 5;
egNumSmall <- nrow(egMatrix);

# debugging code
egSmallMatrix <- egMatrix[1:egNumSmall,];
dataRowNames <- row.names(egSmallMatrix);
resultMatrix <- matrix(nrow=nrow(egSmallMatrix),ncol=0);
#dataColNames <- c("#Feature");
dataColNames <- c();

# Two sample t-test
#! N vs P
# dataTest <- apply(egSmallMatrix, 1, function(x) t.test(x ~ egClass$Class));
dataTest <- apply(egSmallMatrix, 1, function(x) t.test(x ~ egClass$Class[c(indexP, indexN)]));

# retrieved values: t and P-value
dataFTest <- lapply( dataTest, function(x) c(as.numeric(x[1]), as.numeric(x[3])) );
dataFTest <- unlist(dataFTest);
dim(dataFTest) <- c(2, egNumSmall);
dataFTest <- t(dataFTest);

resultMatrix <- cbind(resultMatrix, dataFTest);
dataColNames <- cbind(dataColNames, "Ttest t", "Ttest Pvalue");

# Finalize the result matrix
#resultMatrix <- cbind( dataRowNames, resultMatrix);
#resultMatrix <- rbind( dataColNames, resultMatrix);
colnames(resultMatrix) <- dataColNames;
rownames(resultMatrix) <- dataRowNames;

# heatmap of the top 1000 rows, each statistical test for a column
# We cannot generate the heatmap for complete matrix, due to its memory requirement of 11.1Gb
library(gplots);
egRank <- rank(rank(resultMatrix[,2]));

# heatmap of the top 5 rows, each sample for a column
color.map <- function(tempClass)
{
  if( tempClass=='P' )
    'red'
  else
    'blue';
}
colorCol <- unlist(lapply(egClass$Class, color.map));
egRank <- rank(rank(resultMatrix[,2]));
indexTopRank <- which ( egRank <= 5 ); # top 5 features
resultMatrix[indexTopRank,];


png("tTest_Top5rank_samples_col.png",        
    width = 7*300,        # 7 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size



heatmap.2(as.matrix(egSmallMatrix[indexTopRank,]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);


dev.off()


## histogram

png("tTest_histogram.png",        
    width = 7*300,        # 7 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

hist(resultMatrix[,2], col="gray", labels=TRUE, border="black", breaks=10);

dev.off()
# Names of the top 5 ranked features
egTopFeatures <- names(sort(resultMatrix[indexTopRank, 2]));

heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[1:2],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[2:3],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[3:4],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[4:5],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[1:5],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);

par(mfrow=c(2,2));
plot(as.numeric(egSmallMatrix[egTopFeatures[1],]), as.numeric(egSmallMatrix[egTopFeatures[2],]), col=colorCol);
plot(as.numeric(egSmallMatrix[egTopFeatures[2],]), as.numeric(egSmallMatrix[egTopFeatures[3],]), col=colorCol);
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[4],]), col=colorCol);
plot(as.numeric(egSmallMatrix[egTopFeatures[4],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorCol);
par(mfrow=c(1,1));

# TODO: hclust, kmeans, Density-based, etc

egMatrixGood <- egSmallMatrix[egTopFeatures[3:4],];
egMatrixBad  <- egSmallMatrix[egTopFeatures[4:5],];

par(mfrow=c(1,2));
plot(as.numeric(egMatrixGood[1,]), as.numeric(egMatrixGood[2,]), col=colorCol);
plot(as.numeric(egMatrixBad[1,]), as.numeric(egMatrixBad[2,]), col=colorCol);
par(mfrow=c(1,1));

# Hierarchical clustering
## egMatrixGood
library(ape);
egClusterH <- hclust(dist(t(egMatrixGood)));
parRawMargin <- par("mar");
# bottom, left, top and right margins
par(mar=c(1.2, 2.0, 1.2, 2.0));
par(mfrow=c(1,2));
plot(egClusterH, main="MatrixGood, HCluster");
plot(as.phylo(egClusterH), tip.color=colorCol, main="MatrixGood, HCluster with Color");
par(mfrow=c(1,1));
par(mar=parRawMargin);

egClusterH2 <- cutree(egClusterH, 2);
mapClusterColor <- function(tempClass)
{
    if( tempClass==2 )
        'red'
    else
        'blue';
}
colorCluster <- unlist(lapply(egClusterH2, mapClusterColor));
par(mfrow=c(1,2));
plot(as.numeric(egMatrixGood[1,]), as.numeric(egMatrixGood[2,]), col=colorCol, main="Original annotation" );
plot(as.numeric(egMatrixGood[1,]), as.numeric(egMatrixGood[2,]), col=colorCluster, main="Hierarchical clustering (2)" );
par(mfrow=c(1,1));

egClusterH3 <- cutree(egClusterH, 3);
mapClusterColor <- function(tempClass)
{
    if( tempClass==3 )
        'red'
    else if ( tempClass==2 )
        'green'
    else
        'blue';
}
colorCluster <- unlist(lapply(egClusterH3, mapClusterColor));
par(mfrow=c(1,2));
plot(as.numeric(egMatrixGood[1,]), as.numeric(egMatrixGood[2,]), col=colorCol, main="Original annotation" );
plot(as.numeric(egMatrixGood[1,]), as.numeric(egMatrixGood[2,]), col=colorCluster, main="Hierarchical clustering (3)" );
par(mfrow=c(1,1));



## egMatrixBad
egClusterH <- hclust(dist(t(egMatrixBad)));
parRawMargin <- par("mar");
# bottom, left, top and right margins
par(mar=c(1.2, 2.0, 1.2, 2.0));
par(mfrow=c(1,2));
plot(egClusterH, main="MatrixBad, HCluster");
plot(as.phylo(egClusterH), tip.color=colorCol, main="MatrixBad, HCluster with Color");
par(mfrow=c(1,1));
par(mar=parRawMargin);


# Kmeans
# K-Means Cluster Analysis
tMatrix <- t(egMatrixGood);
tFit <- kmeans(tMatrix, 2) # 2 cluster solution
# get cluster means 
aggregate(tMatrix,by=list(tFit$cluster),FUN=mean)
# append cluster assignment
egMatrixGoodK <- data.frame(tMatrix, tFit$cluster)

mapClusterColor <- function(tempClass)
{
    if( tempClass==1 )
        'red'
    else
        'blue';
}
colorCluster <- unlist(lapply(egMatrixGoodK[,3], mapClusterColor));
par(mfrow=c(1,2));
plot(as.numeric(egMatrixGoodK[,1]), as.numeric(egMatrixGoodK[,2]), col=colorCol, main="Original annotation" );
plot(as.numeric(egMatrixGoodK[,1]), as.numeric(egMatrixGoodK[,2]), col=colorCluster, main="Kmeans clustering (2)" );
par(mfrow=c(1,1));

# 3 clusters
tMatrix <- t(egMatrixGood);
tFit <- kmeans(tMatrix, 3) # 3 cluster solution
# get cluster means 
aggregate(tMatrix,by=list(tFit$cluster),FUN=mean)
# append cluster assignment
egMatrixGoodK <- data.frame(tMatrix, tFit$cluster)

mapClusterColor <- function(tempClass)
{
    if( tempClass==1 )
        'red'
    else if ( tempClass==2 )
        'green'
    else
        'blue';
}
colorCluster <- unlist(lapply(egMatrixGoodK[,3], mapClusterColor));
par(mfrow=c(1,2));
plot(as.numeric(egMatrixGoodK[,1]), as.numeric(egMatrixGoodK[,2]), col=colorCol, main="Original annotation" );
plot(as.numeric(egMatrixGoodK[,1]), as.numeric(egMatrixGoodK[,2]), col=colorCluster, main="Kmeans clustering (3)" );
par(mfrow=c(1,1));


# Model based clustering
# the optimal model according to BIC for EM initialized by hierarchical clustering
# for parameterized Gaussian mixture models.
library(mclust);
tMatrix <- t(egMatrixGood);
# The R function Mclust performs model-based clustering for a range of models
# and a variety of values of k:
tMatrixMCust <- Mclust(tMatrix);
# By default, the models considered are:
#- http://en.wikibooks.org/wiki/Data_Mining_Algorithms_In_R/Clustering/Expectation_Maximization_(EM)
# "EII": spherical, equal volume 
# "VII": spherical, unequal volume 
# "EEI": diagonal, equal volume and shape
# "VEI": diagonal, varying volume, equal shape
# "EVI": diagonal, equal volume, varying shape 
# "VVI": diagonal, varying volume and shape 
# "EEE": ellipsoidal, equal volume, shape, and orientation 
# "EEV": ellipsoidal, equal volume and equal shape
# "VEV": ellipsoidal, equal shape 
# "VVV": ellipsoidal, varying volume, shape, and orientation  

# Plotting the BIC values, and the larger BIC the better
plot(tMatrixMCust, data=tMatrix, what="BIC")
# EII(2) is the best model
# The clustering vector:

tMatrixMCustClass <- tMatrixMCust$classification;
# append cluster assignment
egMatrixGoodK <- data.frame(tMatrix, tMatrixMCustClass)

mapClusterColor <- function(tempClass)
{
    if( tempClass==2 )
        'red'
    else
        'blue';
}
colorCluster <- unlist(lapply(egMatrixGoodK[,3], mapClusterColor));
par(mfrow=c(1,2));
plot(as.numeric(egMatrixGoodK[,1]), as.numeric(egMatrixGoodK[,2]), col=colorCol, main="Original annotation" );
plot(as.numeric(egMatrixGoodK[,1]), as.numeric(egMatrixGoodK[,2]), col=colorCluster, main="Model-based clustering (2)" );
par(mfrow=c(1,1));


# DBSCAN
library(fpc);
tMatrix <- t(egMatrixGood);
tClusterDBSCAN <- dbscan(tMatrix,0.15);

tMatrixMCustClass <- tClusterDBSCAN$cluster;
# append cluster assignment
egMatrixGoodK <- data.frame(tMatrix, tMatrixMCustClass)

mapClusterColor <- function(tempClass)
{
    if( tempClass==0 )
        'red'
    else if ( tempClass==1 )
        'blue'
    else
        'green';
}
colorCluster <- unlist(lapply(egMatrixGoodK[,3], mapClusterColor));
par(mfrow=c(1,2));
plot(as.numeric(egMatrixGoodK[,1]), as.numeric(egMatrixGoodK[,2]), col=colorCol, main="Original annotation" );
plot(as.numeric(egMatrixGoodK[,1]), as.numeric(egMatrixGoodK[,2]), col=colorCluster, main="Model-based clustering (3)" );
par(mfrow=c(1,1));

