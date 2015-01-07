
#***************************************************
#
#
# lian.zou@siat.ac.cn
#
# Visual analysis many kinds of test resuts

#**************************************************

fileMatrix <- read.csv(file="4test-result-GSE10810.csv")
colnum <- ncol(fileMatrix)

resultMatrix <- fileMatrix[,2:colnum]

# 
# # creates a 1500*1500 Images
png("4test-Boxplot.png",        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

# plot
boxplot(resultMatrix[,c(2, 4, 6, 8)]);

dev.off()
