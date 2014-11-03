# *************课堂作业***********************************
#
# 邹炼 ，医工所微创中心，工号：XJY102908
# lian.zou@siat.ac.cn
#
#**********************************************************


#******************* Defined functions*****************************
#lowerTrianglemin
#lowerTrianglemax
#minmaxEuclideanDistance
#******************************************************************



#Dist to Matrix lower triangle of the matrix is used, the rest is ignored
#Get the minimun location in the lower triangle region of the  Matrix,

lowerTrianglemin<-function(matrix){
  dimensions <-dim(matrix)
  rown<- dimensions[1]
  coln<- dimensions[2]
  
  if (rown != coln | rown < 2){
    print("Please ensure the matrix is a square matrix")
    return()
  }
  
  if (rown == 2){
    return(c(1,2))
    
  }
  
  minimun <- matrix[2,1]
  minindex <- c(2,1)
  
  for(i in 3:rown){
  
    jlimit<- i-1
    
    for(j in 1:jlimit){

      if(minimun > matrix[i,j]){
        minimun<- matrix[i,j]
        minindex <- c(i,j)
      }
          
    }
    
    
  }
  
  return(minindex)
  
}



#Dist to Matrix lower triangle of the matrix is used, the rest is ignored
#Get the maximun location in the lower triangle region of the  Matrix,

lowerTrianglemax<-function(matrix){
  dimensions <-dim(matrix)
  rown<- dimensions[1]
  coln<- dimensions[2]
  
  if (rown != coln | rown < 2){
    print("Please ensure the matrix is a square matrix")
    return()
  }
  
  if (rown == 2){
    return(c(1,2))
    
  }
  
  maximun <- matrix[2,1]
  maxindex <- c(2,1)
  
  for(i in 3:rown){
    
    jlimit<- i-1
    
    for(j in 1:jlimit){
      
      if(maximun < matrix[i,j]){
        maximun<- matrix[i,j]
        maxindex <- c(i,j)
      }
      
    }
    
    
  }
  
  return(maxindex)
  
  
}




#Find the Minimun and Maximun in the Distance Matrix
#Return the (minimu,maximu) location:
#As the Euclide Distance, the returned vale point to the sample id.
#Column means the feature as default

minmaxEuclideanDistance <- function(matrix,colfeature=TRUE,print=TRUE) {
  if (colfeature) {
    
    tmatrix <- t(matrix)
    d = dist(tmatrix)
    
    dmatrix <- as.matrix(d)
    
    minindex<-lowerTrianglemin(dmatrix)
    maxindex<-lowerTrianglemax(dmatrix)
    
    minmaxindex<- c(minindex, maxindex)
    return(minmaxindex)
  
  
  } else {
   
    d = dist(matrix)
    
    dmatrix <- as.matrix(d)
    
    
    minindex<-lowerTrianglemin(dmatrix)
    maxindex<-lowerTrianglemax(dmatrix)
    
    minmaxindex<- c(minindex, maxindex)
    return(minmaxindex)
  }
 
}




