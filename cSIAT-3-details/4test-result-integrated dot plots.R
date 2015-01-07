
#****************************************************************
#
# lian.zou@siat.ac.cn
#
# Visual analysis many kinds of test resuts

#**********************************************************
fileMatrix <- read.csv(file="4test-result-GSE10810.csv")
colnum <- ncol(fileMatrix)

resultMatrix <- fileMatrix[,2:colnum]

#integrated dot plots
# # creates a 1500*1500 Images
png("4test-IntegratedDotplots.png",        
    width = 6*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

# integrated dot plots
par(mfrow=c(4,4));
indexCol <- 2*c(1:4);
for ( i in 1:4 )
  for (j in 1:4)
  {
    plot(resultMatrix[,indexCol[i]], resultMatrix[,indexCol[j]])
  }
par(mfrow=c(1,1));


dev.off()
