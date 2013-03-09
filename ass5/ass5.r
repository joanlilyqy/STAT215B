#### Ying Qiao
#### SID: 21412301
#### STAT 215B, Assignment 5, 3/15/2013

rm(list=ls())

#### Simulation
#### Efron (2010) Table 1.2
n <- 10
N <- 1000
# true parameter value
mu <- c(-.81, -.39, -.39, -.08, .69, .71, 1.28, 1.32, 1.89, 4.00)
A.hat <- sum(mu^2) / n
TSE <- (n*A.hat+2)/(A.hat+1)

mu.hat.MLE <- matrix(0, nrow=n, ncol=N)
mu.hat.JS <- matrix(0, nrow=n, ncol=N)

for (i in 1:N){
# random sample
set.seed(i*10)
z <- rnorm(n, mu, 1)
# MLE
mu.hat.MLE[ ,i] <- z
# JS
z.bar <- mean(z)
S <- sum((z - z.bar)^2)
mu.hat.JS[ ,i] <- z.bar + (1 - (n-3)*1/S)*(z - z.bar)
}

MSE.MLE <- rowMeans((mu.hat.MLE - mu)^2)
MSE.MLE
sum(MSE.MLE)
MSE.JS <- rowMeans((mu.hat.JS - mu)^2)
MSE.JS
sum(MSE.JS)




