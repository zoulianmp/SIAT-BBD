
# *************BBE-Final Exam***********************************
#
# Lian Zou, MSI, Study ID£ºXJY102908  
#
# lian.zou@siat.ac.cn
#
# given the samples data
# do many kinds of testï¼? and export the resuts to file

#**********************************************************


# Initialization: loading the dataset
egMatrix <- read.csv("class-matrix-GSE10810.csv", header=TRUE, sep=",", row.names=1);
egClass <- read.csv("class-GSE10810-binary.csv", header=TRUE, sep=",", row.names=1);
indexP <- which(egClass$Class == "P");
indexN <- which(egClass$Class == "N");

#egNumSmall <- 1000; #used for debug 
egNumSmall <- nrow(egMatrix);

egSmallMatrix <- egMatrix[1:egNumSmall,];
dataRowNames <- row.names(egSmallMatrix); #used for resultMatrix rownames
resultMatrix <- matrix(nrow=nrow(egSmallMatrix),ncol=0);

#dataColNames <- c("#Feature");
dataColNames <- c();  #used for resultMatrix colnames






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

# Wilcoxon rank sum test
#! N vs P
# dataTest <- apply(egSmallMatrix, 1, function(x) t.test(x ~ egClass$Class));
dataTest <- apply(egSmallMatrix, 1, function(x) wilcox.test(x ~ egClass$Class[c(indexP, indexN)]));

# retrieved values: t and P-value
dataFTest <- lapply( dataTest, function(x) c(as.numeric(x[1]), as.numeric(x[3])) );
dataFTest <- unlist(dataFTest);
dim(dataFTest) <- c(2, egNumSmall);
dataFTest <- t(dataFTest);

resultMatrix <- cbind(resultMatrix, dataFTest);
dataColNames <- cbind(dataColNames, "Wtest W", "Wtest Pvalue");


# ANOVA test
#! N vs P
# dataTest <- apply(egSmallMatrix, 1, function(x) t.test(x ~ egClass$Class));
### dataTest <- apply(egSmallMatrix, 1, function(x) anovalm((x[indexP], x[indexN])));
dataTest <- apply(egSmallMatrix, 1, function(x) anova(lm(x ~ egClass$Class[c(indexP, indexN)])));

# retrieved values: t and P-value
dataFTest <- lapply( dataTest, function(x) c(as.numeric(unlist(x[4]))[1], as.numeric(unlist(x[5]))[1]) );
dataFTest <- unlist(dataFTest);
dim(dataFTest) <- c(2, egNumSmall);
dataFTest <- t(dataFTest);

resultMatrix <- cbind(resultMatrix, dataFTest);
dataColNames <- cbind(dataColNames, "ANOVA F", "ANOVA Pr(>F)");

# Kruskal-Wallis rank sum test
#! N vs P
# dataTest <- apply(egSmallMatrix, 1, function(x) t.test(x ~ egClass$Class));
dataTest <- apply(egSmallMatrix, 1, function(x) kruskal.test(x ~ egClass$Class[c(indexP, indexN)]));

# retrieved values: t and P-value
dataFTest <- lapply( dataTest, function(x) c(as.numeric(x[1]), as.numeric(x[3])) );
dataFTest <- unlist(dataFTest);
dim(dataFTest) <- c(2, egNumSmall);
dataFTest <- t(dataFTest);

resultMatrix <- cbind(resultMatrix, dataFTest);
dataColNames <- cbind(dataColNames, "KS Chi", "KS Pvalue");

# Finalize the result matrix
#resultMatrix <- cbind( dataRowNames, resultMatrix);
#resultMatrix <- rbind( dataColNames, resultMatrix);
colnames(resultMatrix) <- dataColNames;
rownames(resultMatrix) <- dataRowNames;


write.csv(resultMatrix, "4test-result-GSE10810.csv")
