
# ******************************
#
# 
# lian.zou@siat.ac.cn
#
# Load the t test result from Many Analysis
#
# This Module used for refine tunning the Selected Features
#
#**********************************************************




#*****************************************************
#Initial loading the data, t-test results by zoulian
#******************************************************

egMatrix <- read.csv("class-matrix-GSE10810.csv", header=TRUE, sep=",", row.names=1);

#egNumSmall <- 5;
egNumSmall <- nrow(egMatrix);

# # debugging code
egSmallMatrix <- egMatrix[1:egNumSmall,];


egClass <- read.csv("class-GSE10810-binary.csv", header=TRUE, sep=",", row.names=1);
indexP <- which(egClass$Class == "P");
indexN <- which(egClass$Class == "N");


fileMatrix <- read.csv(file="4test-result-GSE10810.csv")
colnum <- ncol(fileMatrix)

resultMatrix <- fileMatrix[,2:3] #Get the t Test result

dataColNames <- cbind("Ttest t", "Ttest Pvalue");
dataRowNames <- fileMatrix[,1]
  
colnames(resultMatrix) <- dataColNames;
rownames(resultMatrix) <- dataRowNames;






#*****************************************************
#Prepare the Selected Features and FeatureMatrix
#******************************************************

#set all the fine tunned features
RRFSelectedFeatrues<-c('229839_at')

Top2RankFeatures<- c('212653_s_at','202317_s_at')



#features<- t(RRFSelectedFeatrues) #only for RRF selected featureTest

features<-cbind(t(RRFSelectedFeatrues),t(Top2RankFeatures))


featuresMatrix <- matrix(nrow=0,ncol=ncol(egSmallMatrix));

nfeatures <- ncol(features);

for (i in 1:nfeatures){
  
  featureindex<-which(dataRowNames == features[i]);

  featuresMatrix<- rbind(featuresMatrix,egSmallMatrix[featureindex,])
 
}


rownames(featuresMatrix) <- t(features);



#features

#featuresMatrix




#*****************************************************
#Prepare the Analysis Packages and self defined Functions
#******************************************************

# Working
library(MASS)
library(glmnet);#lasso
library(kernlab);#svm
library(rpart);#dtree
library(e1071);#bayes
library(pamr);#pam
library(class); ##K-NN
library(minerva);#mine
#library(FSelector);#best.first
library(RRF);
#library(genefilter);
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


efPerformanceMatrix <- function( kMatrix )
{
  tresultMatrix <- matrix(nrow=0,ncol=5);
  tdataRowNames <- c();
  tdataColNames <- c("Sn", "Sp", "Acc", "Avc", "MCC");
  
  # SVM
  egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "SVM");
  etMeasurements <- egResult$out;
  tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
  tdataRowNames <- c(tdataRowNames, "SVM");
  
  # NBayes
  egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "NBayes");
  etMeasurements <- egResult$out;
  tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
  tdataRowNames <- c(tdataRowNames, "NBayes");
  
  # DTree
  egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "DTree");
  etMeasurements <- egResult$out;
  tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
  tdataRowNames <- c(tdataRowNames, "DTree");
  
  # Lasso
  egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "Lasso");
  etMeasurements <- egResult$out;
  tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
  tdataRowNames <- c(tdataRowNames, "Lasso");
  
  # KNN
  egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "KNN");
  etMeasurements <- egResult$out;
  tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
  tdataRowNames <- c(tdataRowNames, "KNN");
  
  rownames(tresultMatrix) <- tdataRowNames;
  colnames(tresultMatrix) <- tdataColNames;
  
  return ( tresultMatrix );
}



#Generate the Single Feature Performance Matrix
efPerformanceSingleFeature <- function( kMatrix )
{
  
    tresultMatrix <- matrix(nrow=0,ncol=5);
    tdataRowNames <- c();
    tdataColNames <- c("Sn", "Sp", "Acc", "Avc", "MCC");
    
    # SVM
    egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "SVM");
    etMeasurements <- egResult$out;
    tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
    tdataRowNames <- c(tdataRowNames, "SVM");
    
    # NBayes
    egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "NBayes");
    etMeasurements <- egResult$out;
    tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
    tdataRowNames <- c(tdataRowNames, "NBayes");
    
    # DTree
    egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "DTree");
    etMeasurements <- egResult$out;
    tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
    tdataRowNames <- c(tdataRowNames, "DTree");
   
    # KNN
    egResult <- efKFCV(t(kMatrix), egClassLabel, 3, "KNN");
    etMeasurements <- egResult$out;
    tresultMatrix <- rbind(tresultMatrix, etMeasurements[1,]);
    tdataRowNames <- c(tdataRowNames, "KNN");
    
    rownames(tresultMatrix) <- tdataRowNames;
    colnames(tresultMatrix) <- tdataColNames;
    
    return ( tresultMatrix );
    
}



### End of function definition





#**********************************************
#Selected Features Analysis
#**********************************************


color.map <- function(tempClass)
{
  if( tempClass=='P' )
    'red'
  else
    'blue';
}
colorCol <- unlist(lapply(egClass$Class, color.map));



library(gplots);

Selectedfeaturestitle<- paste("RRF Selected 1 + Top rank 2 ----",nfeatures,"Features")

nfeatures


if ( nfeatures > 1 ){
  
    heatmap.2(as.matrix(featuresMatrix), ColSideColors=colorCol, col=redgreen(75), main=Selectedfeaturestitle, margins=c(5,5), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
    set.seed(0);
    egClassLabel <- efClassToInt(as.numeric(egClass$Class));
    resultMatrix <- efPerformanceMatrix(featuresMatrix);
    par(mar=c(3, 2, 1.2, 2)); # bottom, left, top, right
    barplot(resultMatrix, beside=TRUE, col=rainbow(5), main=Selectedfeaturestitle, ylim=c(-0.2, 1.0), xlim=c(0, 36));
    #legend(30, 1.0, rownames(egSmallMatrix[indexTopRank,]), cex=1.0, bty="n", fill=rainbow(5));
    legend(30, 1.0, rownames(resultMatrix), cex=1.0, bty="n", fill=rainbow(5));
    box();
} else if (nfeatures == 1){
    Singlefeaturestitle <- paste("Single Feature Performance : ",features )
    
      
    set.seed(0);
    egClassLabel <- efClassToInt(as.numeric(egClass$Class));
    resultMatrix <- efPerformanceSingleFeature(featuresMatrix);
    par(mar=c(3, 2, 1.2, 2)); # bottom, left, top, right
    barplot(resultMatrix, beside=TRUE, col=rainbow(4), main=Singlefeaturestitle, ylim=c(-0.2, 1.0), xlim=c(0, 36));
    #legend(30, 1.0, rownames(egSmallMatrix[indexTopRank,]), cex=1.0, bty="n", fill=rainbow(5));
    legend(30, 1.0, rownames(resultMatrix), cex=1.0, bty="n", fill=rainbow(5));
    box();
    
    
  
}
  
  















