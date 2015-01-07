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
  }
  return (predy)
}

# SVM
egResult <- efKFCV(t(egSmallMatrix[indexTopRank[c(3,5)],]), egClassLabel, 3, "SVM");
# Sn, Sp, Acc, Avc, MCC
etMeasurements <- egResult$out;


