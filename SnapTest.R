
library(gplots);
par(mfrow=c(3,1));

for ( i in 1:3 )
{
  plot(ecdf(rnorm(10)));
  plot(ecdf(rnorm(10)), add = TRUE, col = 2);
}
plot

par(mfrow=c(1,1));

x3 <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5)
x4 <- c(3, 1, 3, 1, 5, 4, 8, 6, 9, 3, 5)


sumrank<-(rank(x3)+rank(x4))

rank(sumrank)

x1 <- c(3, 1, 4, 15, 1,92,4)
x2 <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5)


(r1 <- rank(x1))

names(x2) <- letters[1:11]
(r2 <- rank(x2)) # ties are averaged

## rank() is "idempotent": rank(rank(x)) == rank(x) :
stopifnot(rank(r1) == r1, rank(r2) == r2)

## ranks without averaging
r.first<-rank(x2, ties.method= "first")  # first occurrence wins
r.random1<-rank(x2, ties.method= "random") # ties broken at random

r.random2<-rank(x2, ties.method= "random") # and again
r.average<-rank(x2, ties.method= "average") # and again


## keep ties ties, no average
(rma <- rank(x2, ties.method= "max"))  # as used classically
(rmi <- rank(x2, ties.method= "min"))  # as in Sports
stopifnot(rma + rmi == round(r2 + r2))
