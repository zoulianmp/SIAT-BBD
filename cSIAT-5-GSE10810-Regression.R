#-------------------------------------------------------------------------------
# R scripts for demonstration of the algorithms taught in the course SIAT-BMI
# Teacher: Fengfeng Zhou
# Email: FengfengZhou@gmail.com
# Web: http://healthinformaticslab.org/course/SIAT-BMI/
# Semester: August 2014
# Version 1.0.1
# Update: 2014-09-01
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

memory.limit(12000);
set.seed(0);

# Initialization: loading the dataset
egMatrix <- read.csv("class-matrix-GSE32175.csv", header=TRUE, sep=",", row.names=1);
egClass <- read.csv("class-GSE32175-binary.csv", header=TRUE, sep=",", row.names=1);
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
indexTopOneKRank <- which ( egRank<=1000 );
indexTopRank <- which ( egRank==1 )[1];
indexRank2 <- which ( egRank==2 )[1];
indexRank3 <- which ( egRank==3 )[1];
indexRank5 <- which ( egRank==5 )[1];
indexRank7 <- which ( egRank==7 )[1];
egValue2 <- as.numeric(egMatrix[indexRank2,]);
egValue3 <- as.numeric(egMatrix[indexRank3,]);
egValue5 <- as.numeric(egMatrix[indexRank5,]);
egValue7 <- as.numeric(egMatrix[indexRank7,]);

# Names of the top 5 ranked features
egTopFeatures <- names(sort(resultMatrix[indexTopRank, 2]));
egTopClass <- as.numeric(egSmallMatrix[indexTopRank,]);
egMatrix <- egSmallMatrix[-indexTopRank,];

indexOrder <- sort(egTopClass, index.return=TRUE);
egTopClass <- egTopClass[indexOrder$ix];
egValue2 <- egValue2[indexOrder$ix];
egValue3 <- egValue3[indexOrder$ix];
egValue5 <- egValue5[indexOrder$ix];
egValue7 <- egValue7[indexOrder$ix];
egColor <- colorCol[indexOrder$ix];

egTopOneKMatrix <- egSmallMatrix[indexTopOneKRank,];



# Linear Regressoin
egFitLinear <- lm(egTopClass~egValue3+egValue5);
egFittedValue <- egFitLinear$fitted.values;
par(mfrow=c(2,2));
plot(egTopClass, col=egColor, main="Original annotations");
plot(egFittedValue, col=egColor, main=paste("Fitted values (3,5) [AIC=", round(AIC(egFitLinear),3), "]", sep=""));
plot(egTopClass, egFittedValue, col=egColor, main="Original vs Fitted  (3,5) values");
egCMatrix <- cbind(egTopClass, egFittedValue);
colnames(egCMatrix) <- c("Raw", "Fitted");
rownames(egCMatrix) <- colnames(egSmallMatrix);
barplot(egCMatrix, beside=TRUE, col=rainbow(5), ylim=c(0.0, 1.2));
legend("topleft", c("Raw", "Fitted"), cex=1.2, bty="n", fill=rainbow(5));
box();
par(mfrow=c(1,1));


# Linear Regressoin
egFitLinear <- lm(egTopClass~egValue3+egValue5+egValue7);
egFittedValue <- egFitLinear$fitted.values;
par(mfrow=c(2,2));
plot(egTopClass, col=egColor, main="Original annotations");
plot(egFittedValue, col=egColor, main=paste("Fitted values (3,5,7) [AIC=", round(AIC(egFitLinear),3), "]", sep=""));
plot(egTopClass, egFittedValue, col=egColor, main="Original vs Fitted  (3,5,7) values");
egCMatrix <- cbind(egTopClass, egFittedValue);
colnames(egCMatrix) <- c("Raw", "Fitted");
rownames(egCMatrix) <- colnames(egSmallMatrix);
barplot(egCMatrix, beside=TRUE, col=rainbow(5), ylim=c(0.0, 1.2));
legend("topleft", c("Raw", "Fitted"), cex=1.2, bty="n", fill=rainbow(5));
box();
par(mfrow=c(1,1));


# Logistic Regressoin
tMin <- min(egTopClass); tMax <- max(egTopClass); etTopClass <- (egTopClass-tMin)/(tMax-tMin);
egFitLinear <- glm(etTopClass~egValue3+egValue5+egValue7, family=binomial("logit"));
egFittedValue <- egFitLinear$fitted.values;
egFittedValue <- egFittedValue*(tMax-tMin)+tMin;
par(mfrow=c(2,2));
plot(egTopClass, col=egColor, main="Original annotations");
plot(egFittedValue, col=egColor, main=paste("Logistic regression values (3,5,7) [AIC=", round(AIC(egFitLinear),3), "]", sep=""));
plot(egTopClass, egFittedValue, col=egColor, main="Original vs Fitted  (3,5,7) values");
egCMatrix <- cbind(egTopClass, egFittedValue);
colnames(egCMatrix) <- c("Raw", "Fitted");
rownames(egCMatrix) <- colnames(egSmallMatrix);
barplot(egCMatrix, beside=TRUE, col=rainbow(5), ylim=c(0.0, 1.2));
legend("topleft", c("Raw", "Fitted"), cex=1.2, bty="n", fill=rainbow(5));
box();
par(mfrow=c(1,1));

# Gaussian Regressoin
egFitLinear <- glm(egTopClass~egValue3+egValue5+egValue7, family=gaussian);
egFittedValue <- egFitLinear$fitted.values;
par(mfrow=c(2,2));
plot(egTopClass, col=egColor, main="Original annotations");
plot(egFittedValue, col=egColor, main=paste("Guassian regression values (3,5,7) [AIC=", round(AIC(egFitLinear),3), "]", sep=""));
plot(egTopClass, egFittedValue, col=egColor, main="Original vs Fitted  (3,5,7) values");
egCMatrix <- cbind(egTopClass, egFittedValue);
colnames(egCMatrix) <- c("Raw", "Fitted");
rownames(egCMatrix) <- colnames(egSmallMatrix);
barplot(egCMatrix, beside=TRUE, col=rainbow(5), ylim=c(0.0, 1.2));
legend("topleft", c("Raw", "Fitted"), cex=1.2, bty="n", fill=rainbow(5));
box();
par(mfrow=c(1,1));



# SELF Regressoin
egSelf <- egTopClass;
egFitLinear <- lm(egTopClass~egSelf);
egFittedValue <- egFitLinear$fitted.values;
par(mfrow=c(2,2));
plot(egTopClass, col=egColor, main="Original annotations");
plot(egFittedValue, col=egColor, main=paste("SELF values [AIC=", round(AIC(egFitLinear),3), "]", sep=""));
plot(egTopClass, egFittedValue, col=egColor, main="Original vs Fitted  (3,5,7) values");
egCMatrix <- cbind(egTopClass, egFittedValue);
colnames(egCMatrix) <- c("Raw", "Fitted");
rownames(egCMatrix) <- colnames(egSmallMatrix);
barplot(egCMatrix, beside=TRUE, col=rainbow(5), ylim=c(0.0, 1.2));
legend("topleft", c("Raw", "Fitted"), cex=1.2, bty="n", fill=rainbow(5));
box();
par(mfrow=c(1,1));


# Regression to the TopOne features
egSelf <- egTopClass;
tSampleNames <- colnames(egTopOneKMatrix);
tFeatureNames <- rownames(egTopOneKMatrix);
tcFeatureNames <- paste("PS", tFeatureNames, sep="");
rownames(egTopOneKMatrix) <- tcFeatureNames;
egListMatrix <- as.data.frame(t(egTopOneKMatrix));
tVariable <- paste("PS", "236431_at", sep="");
# gFitLinear <- lm(PS236431_at~., data=egListMatrix);
egFitLinear <- lm(eval(as.name(tVariable))~., data=egListMatrix);
indexNA <- which ( egFitLinear$coefficients != "NA" );
egCoefficients <- round(egFitLinear$coefficients[indexNA], 5);
egVarNames <- names(egFitLinear$coefficients[indexNA]);
egVarNames <- sub("^PS", "", egVarNames);
egTFormula <- paste( "(", egCoefficients, ")x(", egVarNames, ")", sep="");
egTFormula[1] <- paste("(", egCoefficients[1], ")", sep="");
egNewLineSpace <- 3; # a new line for every two elements
egSpacer <- rep("", length(egTFormula));
for( i in 1:(length(egTFormula)-1) )
{
    if( i %% egNewLineSpace==0 )
    {
        egSpacer[i] <- "\n";
    }
}
egTFormula <- paste( egTFormula, egSpacer, sep="");
egFormula <- paste( egTFormula, collapse=" + ");
egFormula <- paste( "  ", egFormula );
cat(egFormula);

egValueRaw <- as.numeric(egSmallMatrix[indexTopRank,]);
egValueRegression <- as.numeric(1*egCoefficients[1] + apply(egSmallMatrix[egVarNames[-1],]*egCoefficients[-1], 2, sum) );

indexOrderRaw <- sort(egValueRaw, index.return=TRUE);

#plot(egValueRaw, egValueRegression, col=egColor);
plot(egValueRaw[indexOrderRaw$ix], egValueRegression[indexOrderRaw$ix], col=egColor[indexOrderRaw$ix]);
abline(0,1);
text(0.75, 0.95, egFormula, cex=.98);

