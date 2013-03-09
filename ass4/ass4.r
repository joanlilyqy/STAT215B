#### Ying Qiao
#### SID: 21412301
#### STAT 215B, Assignment 4, 3/8/2013

rm(list=ls())
require(MASS)

#### Sqrt of 2-by-2 matrix
sqrt.m22 <- function(M){
  s <- sqrt(det(M))
  t <- sqrt(sum(diag(M)) + 2*s)
  R <- M / t
  diag(R) <- diag(R) + s / t
  R
}


#### estimation parameter definition
beta = 3
sigma2 = 1
rho = 3/4

#### simulation parameter definition
N.MC = 1000
n = 100

#### storing MC results
beta.hat.ols <- rep(0, N.MC)
beta.hat.ivls <- rep(0, N.MC)
sigma2.hat.ivls <- rep(0, N.MC)
sigma2.hat.zz <- rep(0, N.MC)



for (i in 1:N.MC){
#### LM simulation
# random sample generation
set.seed(i)
seeds <- round(runif(3)*1000)

set.seed(seeds[1])
U <- rnorm(n, 0, 1)
set.seed(seeds[2])
V <- rnorm(n, 0, 1)

mu.err <- rep(0, 2)
Sigma.err <- matrix(c(sigma2,rho,rho,sigma2),2,2)
set.seed(seeds[3])
err <- mvrnorm(n, mu.err, Sigma.err)
eps <- err[ ,1]
del <- err[ ,2]

X <- U + 2*V + del
Y <- X*beta + eps

# OLS estimator
ols <- lm(Y ~ X - 1)
#ols$coeff
beta.hat.ols[i] <- ols$coeff

# IVLS (~2SLS) estimator
iv <- lm(X ~ U + V - 1)
X.hat <- iv$fit
ivls <- lm(Y ~ X.hat - 1)
#ivls$coeff
beta.hat.ivls[i] <- ivls$coeff

# IVLS plugin error estimator
SS <- sum((Y - X*ivls$coeff)^2)
#SS / (n - 1)
sigma2.hat.ivls[i] <- SS / (n - 1)

# transformed error estimator
tZ <- t(cbind(U,V))
ZZ.sqinv <- sqrt.m22(solve(tcrossprod(tZ)))
L <- ZZ.sqinv %*% tZ %*% Y
M <- ZZ.sqinv %*% tZ %*% X
trans <- lm(L ~ M - 1)
#sum(trans$res^2)
sigma2.hat.zz[i] <- sum(trans$res^2)

}


#### MC results
## beta
mean(beta.hat.ols)
sd(beta.hat.ols)
sqrt(mean((beta.hat.ols - beta)^2))

mean(beta.hat.ivls)
sd(beta.hat.ivls)
sqrt(mean((beta.hat.ivls - beta)^2))

#hist
xlim=range(beta.hat.ivls);
lwd=2;col="blue"; xlab="beta hat values"; breaks=seq(140,360,10)
par(mfrow=c(2,1),cex.main= 1)
hist(beta.hat.ols, xlim=xlim, col=col,lwd=lwd,main="",xlab=xlab,ylab="OLS estimator")
hist(beta.hat.ivls,xlim=xlim, col=col,lwd=lwd,main="",xlab=xlab,ylab="IVLS estimator")


## sigma2
mean(sigma2.hat.ivls)
sd(sigma2.hat.ivls)

mean(sigma2.hat.zz)
sd(sigma2.hat.zz)

#hist
lwd=2;col="blue"; xlab="sigma2 hat values"; breaks=seq(140,360,10)
par(mfrow=c(2,1),cex.main= 1)
hist(sigma2.hat.ivls, col=col,lwd=lwd,main="",xlab=xlab,ylab="IVLS Plug-in estimator")
hist(sigma2.hat.zz,   col=col,lwd=lwd,main="",xlab=xlab,ylab="Transformed estimator")


