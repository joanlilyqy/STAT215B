#### Ying Qiao
#### SID: 21412301
#### STAT 215B, Assignment 7, 4/19/2013

rm(list=ls())
require(MASS)

#### Test setup
# sizes
p <- 10
n <- 1000
n.cs <- n * 0.8
n.vs <- n - n.cs
# parameters
alpha <- 0
beta <- sample(1:p,p) / 10
del2 <- 0.4^2 #NSR
sig2 <- del2 * sum(beta^2) * n


## simulations
N.rep <- 100
PMSE.ols <- rep(0, N.rep)
PMSE.copas <- rep(0, N.rep)
PMSE.crosval <- rep(0, N.rep)
K.hat <- rep(0, N.rep)
K.tide <- rep(0, N.rep)

sig2.hat <- rep(0, N.rep)
Vnorm <- rep(0, N.rep)
covXnorm <- rep(0, N.rep)


for (i in 1:N.rep){
#### Sample
# linear model
my.model <- function(x) alpha+crossprod(beta,x)
# CS
X.cs <- matrix(rnorm(n.cs*p), nrow=n.cs, ncol=p)
X.cs <- t(t(X.cs)-colMeans(X.cs)) 
y.cs <- rnorm(n.cs, mean = apply(X.cs, 1, my.model), sd = sqrt(sig2))
V <- 1/n.cs * crossprod(X.cs)
Vnorm[i] <- norm(V, 'F')
covXnorm[i] <- norm(cov(X.cs), 'F')
# VS
X.vs <- mvrnorm(n.vs, rep(0,p), V) #copas
#X.vs <- mvrnorm(n.vs, rep(0,p), V+diag(V)*0.2) #non-copas
X.vs <- t(t(X.vs)-colMeans(X.vs))
y.vs <- rnorm(n.vs, mean = apply(X.vs, 1, my.model), sd = sqrt(sig2))


#### Multiple linear regression on CS
# OLS
ols <- lm(y.cs ~ X.cs)
alpha.ols <- ols$coeff[1]
beta.ols <- ols$coeff[-1]

# K.hat (Copas 1983)
nu <- n.cs - (p+1)
sig2.hat[i] <- sum((ols$res)^2)/nu
bVb <- t(beta.ols) %*% V %*% beta.ols
K.hat[i] <- 1 - (p-2)*sig2.hat[i]*nu/(n.cs*(nu+2)*bVb)

# K.tide (cross validation, ten-fold)
K.cv <- rep(0, 10)
ids <- 1:n.cs
for (j in 1:10){
  # ten-fold
  tt <- sample(ids, n.cs/10)
  Xcv <- X.cs[-tt, ]
  ycv <- y.cs[-tt]
  Xtt <- X.cs[tt, ]
  ytt <- y.cs[tt]
  ids <- ids[-tt]

  # fit
  mm <- lm(ycv ~ Xcv)
  a <- mm$coeff[1]
  b <- mm$coeff[-1]
  y.fit <- apply(Xtt, 1, function(x){a+crossprod(b,x)})
  K.cv[j] <- cov(ytt,y.fit)/var(y.fit)
}
K.tide[i] <- mean(abs(K.cv))


# PMSE (VS)
pred.ols <- function(x) alpha.ols + crossprod(beta.ols, x)
pred.copas <- function(x) alpha.ols + crossprod(K.hat[i]*beta.ols, x)
pred.crosval <- function(x) alpha.ols + crossprod(K.tide[i]*beta.ols, x)
y.ols <- apply(X.vs, 1, pred.ols)
y.copas <- apply(X.vs, 1, pred.copas)
y.crosval <- apply(X.vs, 1, pred.crosval)
PMSE.ols[i] <- mean((y.vs - y.ols)^2)
PMSE.copas[i] <- mean((y.vs - y.copas)^2)
PMSE.crosval[i] <- mean((y.vs - y.crosval)^2)
}

mean(K.hat)
mean(K.tide)

t.test(PMSE.copas, PMSE.ols, alter = 'less')
t.test(PMSE.copas, PMSE.crosval, alter = 'less')
#t.test(PMSE.copas, PMSE.crosval, alter = 'greater') #non-copas

######################
#check
######################
# randomness
hist(X.cs)

# sig2
t.test(sig2.hat, mu=sig2)

# V vs. cov
t.test(Vnorm, covXnorm)
