


x <- matrix(rnorm(12), nrow=2)


## plot
par(mfrow=c(1,2))

plot(x[1,], x[2,],col.lab = "thistle")

#title(main = "The Minimum of Euclidean Distaces" ,  xlab = "GGAGOG", ylab = "SGGSGGGG")

abline(0,1)

plot(x[1,], x[2,],col.lab = "thistle")

title(main = "The Maximum of Euclidean Distaces" ,  xlab = "GGAGOG", ylab = "SGGSGGGG")

abline(0,1)

par(mfrow=c(1,1))
