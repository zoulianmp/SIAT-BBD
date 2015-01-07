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
indexTopRank <- which ( egRank <= 50 ); # top 50 features
resultMatrix[indexTopRank,];
heatmap.2(as.matrix(egSmallMatrix[indexTopRank,]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);

## histogram
hist(resultMatrix[,2], col="gray", labels=TRUE, border="black", breaks=10);

# Names of the top 5 ranked features
egTopFeatures <- names(sort(resultMatrix[indexTopRank, 2]));

heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[1:2],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[2:3],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[3:4],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[4:5],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[c(1,5)],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);

par(mfrow=c(2,2));
plot(as.numeric(egSmallMatrix[egTopFeatures[1],]), as.numeric(egSmallMatrix[egTopFeatures[2],]), col=colorCol);
plot(as.numeric(egSmallMatrix[egTopFeatures[2],]), as.numeric(egSmallMatrix[egTopFeatures[3],]), col=colorCol);
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[4],]), col=colorCol);
plot(as.numeric(egSmallMatrix[egTopFeatures[4],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorCol);
par(mfrow=c(1,1));

par(mfrow=c(1,1));
heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[c(3,5)],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorCol);

# Working
library(MASS)
library(glmnet);#lasso
library(kernlab);#svm
library(rpart);#dtree
library(e1071);#bayes
library(pamr);#pam
library(class); ##K-NN
library(minerva);#mine
library(FSelector);#best.first
library(RRF);
library(genefilter);
library(caret); # createFolds()

### Definitions of all the functions
#! k-Fold Cross Validation
efKFCV <- function(xx,yy,nfold,method)
{
    num_tp=num_fp=num_fn=num_tn=0
    index=NULL
    predy=NULL
    id=createFolds(yy, k = nfold, list = TRUE, returnTrain = T)
    rawdata=cbind(xx,yy)
    n=nrow(rawdata)
    p=ncol(rawdata)
    
    tPrediction <- rep(-1, nrow(rawdata));
    
    for (i in 1:nfold){
        # print(paste("Fold",i,sep=' '))
        index=id[[i]]
        y_train=rawdata[index,p]
        y_test=rawdata[-index,p]
        x_train=matrix(rawdata[index,-p],nrow=length(index))
        x_test=matrix(rawdata[-index,-p],nrow=(n-length(index)))
        
        predy <- efClassifier(x_train,y_train,x_test,method)
        
        num_tp[i]=sum(y_test==1 & predy==1)
        num_fn[i]=sum(y_test==1 & predy==0)
        num_fp[i]=sum(y_test==0 & predy==1)
        num_tn[i]=sum(y_test==0 & predy==0)
        
        tPrediction[-index] <- predy;
    }
    se=sum(num_tp)/sum(yy==1)
    sp=sum(num_tn)/sum(yy==0)
    acc=sum(num_tp+num_tn)/length(yy)
    avc=(se+sp)/2
    mcc=(sum(num_tp)*sum(num_tn)-sum(num_fp)* sum(num_fn))/
        (sqrt((sum(num_tp)+sum(num_fp))*(sum(num_tp)+sum(num_fn))*(sum(num_tn)+sum(num_fp))*(sum(num_tn)+sum(num_fn))))
    out=round(cbind(se,sp,acc,avc,mcc),3)
    return(list(out=out, prediction=tPrediction))
}
#! Standard IO of all the investigated classifiers
efClassifier <- function(x_train,y_train,x_test,method)
{
    if (method=="SVM"){
        fit=ksvm(x_train,y_train,type="C-svc",kernel="rbfdot")
        predy=predict(fit,x_test)
    } else { 
        if (method=="NBayes"){
            colnames(x_train)=NULL
            colnames(x_test)=NULL
            data_train=data.frame(ex=x_train,ey=as.factor(y_train))
            data_test=with(data_train,data.frame(ex=x_test))
            fit <- naiveBayes(ey~.,data=data_train)
            predy=predict(fit, data_test, type="class")
        } else {
            if (method=="DTree"){
                colnames(x_train)=NULL
                colnames(x_test)=NULL
                data_train=data.frame(ex=x_train,ey=as.factor(y_train))
                data_test=with(data_train,data.frame(ex=x_test))
                fit <- rpart(ey~.,data=data_train)
                predy=predict(fit, data_test, type="class")
            } else {
                if (method=="Lasso"){
                    cv.fit <- cv.glmnet(x_train, y_train, family = "binomial")
                    fit <- glmnet(x_train, y_train, family = "binomial")
                    pfit = predict(fit,x_test,s = cv.fit$lambda.min,type="response")
                    predy<-ifelse(pfit>0.5,1,0)
                } else { 
                    if (method=="KNN"){
                        predy<-knn1(x_train,x_test,y_train)
                    }
                }    
            }
        }
    }
    return (predy)
}

efClassToInt<-function(classes)
{
    levelsClass <- sort(levels(as.factor(classes)));
    #for(i in levelsClass)
    for(i in levelsClass)
    {
        classes<-replace(classes,classes==i,match(i,levelsClass)-1);
        #classes <- replace(classes, classes==levelsClass[i], i-1);
    }
    classes <- (as.numeric(classes));
    return (classes);
}

efBinaryPerformance <- function(tClass, tPrediction)
{
    tTP <- sum(tClass==1 & tPrediction==1 );
    tFN <- sum(tClass==1 & tPrediction==0 );
    tFP <- sum(tClass==0 & tPrediction==1 );
    tTN <- sum(tClass==0 & tPrediction==0 );
    tSn <- tTP/(tTP+tFN);
    tSp <- tTN/(tTN+tFP);
    tAcc <- (tTP+tTN)/(tTP+tFN+tTN+tFP);
    tAvc <- (tSn+tSp)/2;
    tMCC <- (tTP*tTN-tFP*tFN)/sqrt((tTP+tFP)*(tTP+tFN)*(tTN+tFP)*(tTN+tFN));
    return(round(cbind(tSn, tSp, tAcc, tAvc, tMCC), 3));
}
### End of function definition


# Main classification functions

#memory.limit(4000)
set.seed(1);

# if egClass$Class is a "factor", use "as.numeric(egClass$Class)" 
egClassLabel <- efClassToInt(as.numeric(egClass$Class));
## egMap <- c("N", "P");

heatmap.2(as.matrix(egSmallMatrix[egTopFeatures[c(3,5)],]), ColSideColors=colorCol, col=redgreen(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);

resultMatrix <- matrix(nrow=0,ncol=5);
dataRowNames <- c();
dataColNames <- c("Sn", "Sp", "Acc", "Avc", "MCC");

# SVM
egResult <- efKFCV(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, 3, "SVM");
etMeasurements <- egResult$out;
eqPrediction   <- egResult$prediction;
#etPrediction <- efClassToInt(eqPrediction);
etPrediction <- eqPrediction;

egResult <- efClassifier(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, t(egSmallMatrix[indexTopRank[c(3,5)],]), "SVM");
qPrediction   <- egResult;
#tPrediction <- efClassToInt(qPrediction);
tPrediction <- qPrediction;
tMeasurements <- efBinaryPerformance(egClassLabel, tPrediction);

# dotplot
par(mfrow=c(1,3));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorCol, main="Original annotation");
color.map <- function(tempClass)
{
    if( tempClass==1 )
        'red'
    else
        'blue';
}
colorPrediction <- unlist(lapply(etPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorPrediction, main=paste("Prediction [SVM] [3FCV-Acc=", etMeasurements[1,"acc"], "]", sep=""));
colorSelfPrediction <- unlist(lapply(tPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorSelfPrediction, main=paste("Prediction [SVM] [Self-Acc=", tMeasurements[1,"tAcc"], "]", sep=""));
par(mfrow=c(1,1));

resultMatrix <- rbind(resultMatrix, etMeasurements[1,]);
dataRowNames <- rbind(dataRowNames, "SVM");


# NBayes
egResult <- efKFCV(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, 3, "NBayes");
etMeasurements <- egResult$out;
eqPrediction   <- egResult$prediction;
etPrediction <- efClassToInt(eqPrediction);

egResult <- efClassifier(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, t(egSmallMatrix[indexTopRank[c(3,5)],]), "NBayes");
qPrediction   <- egResult;
#tPrediction <- efClassToInt(qPrediction);
tPrediction <- qPrediction;
tMeasurements <- efBinaryPerformance(egClassLabel, tPrediction);

# dotplot
par(mfrow=c(1,3));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorCol, main="Original annotation");
color.map <- function(tempClass)
{
    if( tempClass==1 )
        'red'
    else
        'blue';
}
colorPrediction <- unlist(lapply(etPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorPrediction, main=paste("Prediction [NBayes] [3FCV-Acc=", etMeasurements[1,"acc"], "]", sep=""));
colorSelfPrediction <- unlist(lapply(tPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorSelfPrediction, main=paste("Prediction [NBayes] [Self-Acc=", tMeasurements[1,"tAcc"], "]", sep=""));
par(mfrow=c(1,1));

resultMatrix <- rbind(resultMatrix, etMeasurements[1,]);
dataRowNames <- rbind(dataRowNames, "NBayes");



# DTree
egResult <- efKFCV(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, 3, "DTree");
etMeasurements <- egResult$out;
eqPrediction   <- egResult$prediction;
etPrediction <- efClassToInt(eqPrediction);

egResult <- efClassifier(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, t(egSmallMatrix[indexTopRank[c(3,5)],]), "DTree");
qPrediction   <- egResult;
#tPrediction <- efClassToInt(qPrediction);
tPrediction <- qPrediction;
tMeasurements <- efBinaryPerformance(egClassLabel, tPrediction);

# dotplot
par(mfrow=c(1,3));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorCol, main="Original annotation");
color.map <- function(tempClass)
{
    if( tempClass==1 )
        'red'
    else
        'blue';
}
colorPrediction <- unlist(lapply(etPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorPrediction, main=paste("Prediction [DTree] [3FCV-Acc=", etMeasurements[1,"acc"], "]", sep=""));
colorSelfPrediction <- unlist(lapply(tPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorSelfPrediction, main=paste("Prediction [DTree] [Self-Acc=", tMeasurements[1,"tAcc"], "]", sep=""));
par(mfrow=c(1,1));

resultMatrix <- rbind(resultMatrix, etMeasurements[1,]);
dataRowNames <- rbind(dataRowNames, "DTree");



# Lasso
egResult <- efKFCV(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, 3, "Lasso");
etMeasurements <- egResult$out;
eqPrediction   <- egResult$prediction;
etPrediction <- efClassToInt(eqPrediction);

egResult <- efClassifier(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, t(egSmallMatrix[indexTopRank[c(3,5)],]), "Lasso");
qPrediction   <- egResult;
#tPrediction <- efClassToInt(qPrediction);
tPrediction <- qPrediction;
tMeasurements <- efBinaryPerformance(egClassLabel, tPrediction);

# dotplot
par(mfrow=c(1,3));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorCol, main="Original annotation");
color.map <- function(tempClass)
{
    if( tempClass==1 )
        'red'
    else
        'blue';
}
colorPrediction <- unlist(lapply(etPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorPrediction, main=paste("Prediction [Lasso] [3FCV-Acc=", etMeasurements[1,"acc"], "]", sep=""));
colorSelfPrediction <- unlist(lapply(tPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorSelfPrediction, main=paste("Prediction [Lasso] [Self-Acc=", tMeasurements[1,"tAcc"], "]", sep=""));
par(mfrow=c(1,1));

resultMatrix <- rbind(resultMatrix, etMeasurements[1,]);
dataRowNames <- rbind(dataRowNames, "Lasso");



# KNN
egResult <- efKFCV(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, 3, "KNN");
etMeasurements <- egResult$out;
eqPrediction   <- egResult$prediction;
etPrediction <- efClassToInt(eqPrediction);

egResult <- efClassifier(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, t(egSmallMatrix[indexTopRank[c(3,5)],]), "KNN");
qPrediction   <- egResult;
#tPrediction <- efClassToInt(qPrediction);
tPrediction <- qPrediction;
tMeasurements <- efBinaryPerformance(egClassLabel, tPrediction);

# dotplot
par(mfrow=c(1,3));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorCol, main="Original annotation");
color.map <- function(tempClass)
{
    if( tempClass==1 )
        'red'
    else
        'blue';
}
colorPrediction <- unlist(lapply(etPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorPrediction, main=paste("Prediction [KNN] [3FCV-Acc=", etMeasurements[1,"acc"], "]", sep=""));
colorSelfPrediction <- unlist(lapply(tPrediction, color.map));
plot(as.numeric(egSmallMatrix[egTopFeatures[3],]), as.numeric(egSmallMatrix[egTopFeatures[5],]), col=colorSelfPrediction, main=paste("Prediction [KNN] [Self-Acc=", tMeasurements[1,"tAcc"], "]", sep=""));
par(mfrow=c(1,1));

resultMatrix <- rbind(resultMatrix, etMeasurements[1,]);
dataRowNames <- rbind(dataRowNames, "KNN");

rownames(resultMatrix) <- dataRowNames;
colnames(resultMatrix) <- dataColNames;

barplot(resultMatrix, beside=TRUE, col=rainbow(5), ylim=c(-0.2, 1.0));
legend("topright", dataRowNames, cex=1.2, bty="n", fill=rainbow(5));
box();

vRange <- range(resultMatrix);
plot(resultMatrix["SVM",], type="o", col="blue", ylim=vRange, axes=FALSE, ann=FALSE);
axis(1, at=1:5, lab=dataColNames);
axis(2, las=1, at=0:vRange[2]);
box();
lines(resultMatrix["NBayes",], type="o", pch=22, lty=1, col="red");
lines(resultMatrix["DTree",], type="o", pch=15, lty=1, col="purple");
lines(resultMatrix["Lasso",], type="o", pch=25, lty=1, col="green");
lines(resultMatrix["KNN",], type="o", pch=18, lty=1, col="black");
title(main="Comparison of 3-fold cross validation of the five algorithms", col.main="red", font.main=4);
title(xlab="Value", col.lab=rgb(0,0.5,0));
title(ylab="Algorithm", col.lab=rgb(0,0.5,0));
legend(2, 0.4, dataRowNames, cex=1.2, col=c("blue", "red", "purple", "green", "black"), pch=21:22, lty=1);



