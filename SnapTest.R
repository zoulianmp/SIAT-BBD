source("./zoulianfuncs/mysummary.R")


x <- matrix(c(1, 2, 0, 4), nrow = 2, ncol = 2)

which(x == min(x), arr.ind = TRUE)



x <- c(1:4, 0:5, 11)
which.min(x)
which.max(x)


which()

# invoking the function 
set.seed(1234)
x <- rpois(500, 4) 
y <- mysummary(x)
Median= 4
MAD= 1.4826 
# y$center is the median (4) 
# y$spread is the median absolute deviation (1.4826)

y <- mysummary(x, npar=FALSE, print=FALSE)
# no output 
# y$center is the mean (4.052)
# y$spread is the standard deviation (2.01927)



x <- matrix(rnorm(12), nrow=3)
d1 = dist(x)
d2 = dist(x, diag = TRUE)
d3 = dist(x, upper = TRUE)


which.min()
which.max(x)

m <- as.matrix(dist(x))
d <- as.dist(m)
stopifnot(d == dist(x))
#names(d) <- LETTERS[1:5]
print(d, digits = 3)
