
# ******************************
#
# 
# lian.zou@siat.ac.cn
#
# Visual analysis many kinds of test resuts

#**********************************************************

fileMatrix <- read.csv(file="4test-result-GSE10810.csv")
colnum <- ncol(fileMatrix)

resultMatrix <- fileMatrix[,2:colnum]


## histogram
# # creates a 1500*1500 Images
png("4test-histgram.png",        
    width = 6*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size


par(mfrow=c(2,2));
hist(resultMatrix[,2], col="gray", labels=TRUE, border="black", breaks=10, main="Two sample t-test", margins= c(10,10));
hist(resultMatrix[,4], col="gray", labels=TRUE, border="black", breaks=10, main="Wilcoxon rank sum test", margins= c(10,10));
hist(resultMatrix[,6], col="gray", labels=TRUE, border="black", breaks=10, main="ANOVA test", margins= c(10,10));
hist(resultMatrix[,8], col="gray", labels=TRUE, border="black", breaks=10, main="Kruskal-Wallis rank sum test", margins= c(10,10));
par(mfrow=c(1,1));

dev.off()


