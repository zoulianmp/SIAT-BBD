
# ******************************
#
# 
# lian.zou@siat.ac.cn
#
# Load the t test result for Many Analysis

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




#*****************************************************
#Top x Features heatmap 
#******************************************************


dataFTest <-resultMatrix #assign the t test result

#******************
#***

egRankingFeatureNumber <- 10; #The Top X value
egMaxFeatureNumber <- 50;

#***
#*******************



color.map <- function(tempClass)
{
  if( tempClass=='P' )
    'red'
  else
    'blue';
}
colorCol <- unlist(lapply(egClass$Class, color.map));



library(gplots);

topfeaturestitle<- paste("T-test Top",egRankingFeatureNumber,"Features")

egRank <- rank(dataFTest[,2]);
indexTopRank <- which ( egRank <= egRankingFeatureNumber ); # top 10 features
heatmap.2(as.matrix(egSmallMatrix[indexTopRank,]), ColSideColors=colorCol, col=redgreen(75), main=topfeaturestitle, margins=c(10,7), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
set.seed(0);
egClassLabel <- efClassToInt(as.numeric(egClass$Class));
resultMatrix <- efPerformanceMatrix(egSmallMatrix[indexTopRank,]);
par(mar=c(3, 2, 1.2, 2)); # bottom, left, top, right
barplot(resultMatrix, beside=TRUE, col=rainbow(5), main=topfeaturestitle, ylim=c(-0.2, 1.0), xlim=c(0, 36));
#legend(30, 1.0, rownames(egSmallMatrix[indexTopRank,]), cex=1.0, bty="n", fill=rainbow(5));
legend(30, 1.0, rownames(resultMatrix), cex=1.0, bty="n", fill=rainbow(5));
box();




#*****************************************************
#Error and try to find the fewer Features for SVM and NBayes 
#******************************************************

#Generate the SVM and NBayes Performance Matrix
efPerformanceSVMNBayes <- function( kMatrix )
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
  
  
  rownames(tresultMatrix) <- tdataRowNames;
  colnames(tresultMatrix) <- tdataColNames;
  
  return ( tresultMatrix );
}



#Define the function to get sn,sp,mcc of SVM and NBayes
#Use for t Test, get the top nfeatures Performance
#give the For 1:nfeatures loop results
tTestSVMNBayesPerformanceMatrix <- function(ttestout,matrix,nfeatures=10)
{
  tresultMatrix <- matrix(nrow=0,ncol=6);
  tdataRowNames <- c();
  tdataColNames <- c("SVM-sn","SVM-sp","SVM-mcc","NBayes-sn","NBayes-sp","NBayes-mcc");
  
  egRank <- rank(ttestout[,2]);

   
  for (i in 1:nfeatures){
     indexTopRank <- which ( egRank <= i ); # top 10 features
     resultMatrix <- efPerformanceSVMNBayes(matrix[indexTopRank,]);
     svmnbayesperform <-cbind(t(resultMatrix[1,-(3:4)]),t(resultMatrix[2,-(3:4)]))
  
     tresultMatrix <- rbind(tresultMatrix, svmnbayesperform);
         
     tdataRowNames <- c(tdataRowNames, paste(i,"Features"));
  }

  rownames(tresultMatrix) <- tdataRowNames;
  colnames(tresultMatrix) <- tdataColNames;

  return ( tresultMatrix );
  
}
#End of function define

#Do the finding job

out<-tTestSVMNBayesPerformanceMatrix(dataFTest,egSmallMatrix,egRankingFeatureNumber) 

# Plot the performance curves


par(mfrow=c(3,2));
x<-c(1:10)

svmindex<-c(1:3);
nbaysindex<-3+c(1:3);

label<-colnames(out)

for ( i in svmindex )
{
  j<-nbaysindex[i];
  
  plot(out[,i],type="b",main=paste(label[i],"vs Features Number"),col = 2);
  plot(out[,j], type="b", main=paste(label[j],"vs Features Number"),col = 4);
}

par(mfrow=c(1,1));




#*****************************************************
#3 Folder cross validation
#******************************************************

# Group optimization: PAM
#! N vs P


set.seed(0);
PAMdata <- list(x=as.matrix(egSmallMatrix), genenames=paste("g",dataRowNames), geneid=paste("g",dataRowNames), y=factor(egClass$Class));
PAMtrain<- pamr.train(PAMdata);
PAMcv <- pamr.cv(PAMtrain,PAMdata);

PAMpdex1<-max(which(PAMcv$error==min(PAMcv$error)));
PAMpdex2<-max(which(PAMcv$size!=0));
PAMpdex<-ifelse (PAMcv$size[PAMpdex1]==0,PAMpdex2,PAMpdex1);
PAMcvthreshold<-PAMcv$threshold[PAMpdex];

PAMglist <- pamr.listgenes(PAMtrain, PAMdata, PAMcvthreshold, fitcv=PAMcv);
egChosenFeatures <- sub("^g ", "", PAMglist[,1]);
if( length(egChosenFeatures)>egMaxFeatureNumber )
{
  egChosenFeatures <- egChosenFeatures[1:egMaxFeatureNumber];
}


heatmap.2(as.matrix(egSmallMatrix[egChosenFeatures,]), ColSideColors=colorCol, col=redgreen(75), main=paste("PAM chosen ", length(egChosenFeatures), " features"), margins=c(10,7), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
set.seed(0);
egClassLabel <- efClassToInt(as.numeric(egClass$Class));
resultMatrix <- efPerformanceMatrix(egSmallMatrix[egChosenFeatures,]);
par(mar=c(3, 2, 1.2, 2)); # bottom, left, top, right
barplot(resultMatrix, beside=TRUE, col=rainbow(5), main=paste("PAM chosen ", length(egChosenFeatures), " features"), ylim=c(-0.2, 1.0), xlim=c(0, 36));
legend(30, 1.0, rownames(resultMatrix), cex=1.0, bty="n", fill=rainbow(5));
box();



# Group optimization: RRF
#! N vs P


set.seed(0);
RRFrrf <- RRF(as.matrix(t(egSmallMatrix)),y=as.factor(egClass$Class));
RRFimp<-RRFrrf$importance;
RRFimp<-RRFimp[,"MeanDecreaseGini"];
fRRF <- which(RRFimp>0);  ##FS index
egChosenFeatures <- dataRowNames[fRRF];
if( length(egChosenFeatures)>egMaxFeatureNumber )
{
  egChosenFeatures <- egChosenFeatures[1:egMaxFeatureNumber];
}


if (length(egChosenFeatures) > 1 ){
  
   heatmap.2(as.matrix(egSmallMatrix[egChosenFeatures,]), ColSideColors=colorCol, col=redgreen(75), main=paste("RRF chosen ", length(egChosenFeatures), " features"), margins=c(10,7), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
   set.seed(0);
   egClassLabel <- efClassToInt(as.numeric(egClass$Class));
   resultMatrix <- efPerformanceMatrix(egSmallMatrix[egChosenFeatures,]);
   par(mar=c(3, 2, 1.2, 2)); # bottom, left, top, right
   barplot(resultMatrix, beside=TRUE, col=rainbow(5), main=paste("RRF chosen ", length(egChosenFeatures), " features"), ylim=c(-0.2, 1.0), xlim=c(0, 36));
   legend(30, 1.0, rownames(resultMatrix), cex=1.0, bty="n", fill=rainbow(5));
   box();

}else if (length(egChosenFeatures)== 1){
    set.seed(0);
    egClassLabel <- efClassToInt(as.numeric(egClass$Class));
    resultMatrix <- efPerformanceSingleFeature(egSmallMatrix[egChosenFeatures,]);
    par(mar=c(3, 2, 1.2, 2)); # bottom, left, top, right
    barplot(resultMatrix, beside=TRUE, col=rainbow(4), main=paste("RRF chosen ", length(egChosenFeatures), " features"), ylim=c(-0.2, 1.0), xlim=c(0, 36));
    legend(30, 1.0, rownames(resultMatrix), cex=1.0, bty="n", fill=rainbow(5));
    box();
}







# Group optimization: Lasso
#! N vs P
set.seed(0);

#storage.mode£¨as.matrix(t(egSmallMatrix))£©

cv.fit <- cv.glmnet(as.matrix(t(egSmallMatrix)), y=as.factor(egClass$Class), family = "binomial")
lambda<- cv.fit$lambda.min
fit <- glmnet(as.matrix(t(egSmallMatrix)), y=as.factor(egClass$Class),, family = "binomial")
Coef <- coef(fit, s = lambda)[-1]
Active.Index <- which(Coef!= 0)
Active.Coef <- Coef[Active.Index]   

egChosenFeatures <- dataRowNames[Active.Index];
if( length(egChosenFeatures)>egMaxFeatureNumber )
{
  egChosenFeatures <- egChosenFeatures[1:egMaxFeatureNumber];
}

heatmap.2(as.matrix(egSmallMatrix[egChosenFeatures,]), ColSideColors=colorCol, col=redgreen(75), main=paste("Lasso chosen ", length(egChosenFeatures), " features"), margins=c(10,7), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1);
set.seed(0);
egClassLabel <- efClassToInt(as.numeric(egClass$Class));
resultMatrix <- efPerformanceMatrix(egSmallMatrix[egChosenFeatures,]);
par(mar=c(3, 2, 1.2, 2)); # bottom, left, top, right
barplot(resultMatrix, beside=TRUE, col=rainbow(5), main=paste("Lasso chosen ", length(egChosenFeatures), " features"), ylim=c(-0.2, 1.0), xlim=c(0, 36));
legend(30, 1.0, rownames(resultMatrix), cex=1.0, bty="n", fill=rainbow(5));
box();


