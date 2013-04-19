#### Ying Qiao
#### SID: 21412301
#### STAT 215B, Assignment 7, 4/19/2013

rm(list=ls())
require(MASS)

#### Test setup
# sizes
p <- 50
n <- 500
n.cs <- n * 0.8
n.vs <- n - n.cs
# parameters
alpha <- 0.1
beta <- sample(1:p,p) / 10
del2 <- 0.1^2 #NSR
sig2 <- del2 * sum(beta^2) * n.cs
# model
my.model <- function(x, a,b,K) a+K*crossprod(b,x)
  
## simulations
N.rep <- 1000
PMSE<-list()
PMSE$ols <- rep(0, N.rep)
PMSE$copas <- rep(0, N.rep)
PMSE$crosval <- rep(0, N.rep)
K<-list()
K$copas <- rep(0, N.rep)
K$crosval <- rep(0, N.rep)

sig2hat <- rep(0, N.rep)
Vnorm <- rep(0, N.rep)
covXnorm <- rep(0, N.rep)


for (i in 1:N.rep){
#### Sample
# cov of X
r <- runif(1)
vi <- rep(1:p, each=p)
vj <- rep(1:p, p)
cov <- matrix(sapply(1:(p*p), function(i){r^(abs(vi[i]-vj[i]))}), nrow=p)
# CS
X.cs <- mvrnorm(n.cs, rep(0,p), cov)
X.cs <- t(t(X.cs)-colMeans(X.cs))
V <- 1/n.cs * crossprod(X.cs)
Vnorm[i] <- norm(V, 'F')
covXnorm[i] <- norm(cov(X.cs), 'F')

eps <- rnorm(n, mean=0, sd=sqrt(sig2)) #copas
# non-copas
#eps <- rnorm(n/2, mean=-3*sqrt(sig2), sd=sqrt(sig2))
#eps <- c(eps, rnorm(n/2, mean=3*sqrt(sig2), sd=sqrt(sig2)))

ii.cs <- sample(1:n, n.cs)
y.cs <- apply(X.cs, 1, my.model, a=alpha,b=beta,K=1) + eps[ii.cs]
# VS
X.vs <- mvrnorm(n.vs, rep(0,p), cov)
X.vs <- t(t(X.vs)-colMeans(X.vs))
y.vs <- apply(X.vs, 1, my.model, a=alpha,b=beta,K=1) + eps[-ii.cs]


#### Multiple linear regression on CS
# OLS
ols <- lm(y.cs ~ X.cs)
alpha.ols <- ols$coeff[1]
beta.ols <- ols$coeff[-1]

# K.hat (Copas 1983)
nu <- n.cs - (p+1)
sig2hat[i] <- sum((ols$res)^2)/nu
bVb <- t(beta.ols) %*% V %*% beta.ols
K$copas[i] <- 1 - (p-2)*sig2hat[i]*nu/(n.cs*(nu+2)*bVb)

# K.tide (cross validation, k-fold)
k=10
K.cv <- rep(0, k)
pmse.cv <- rep(0, k)
ids <- 1:n.cs
for (j in 1:k){
  # ten-fold
  tt <- sample(ids, n.cs/k)
  Xcv <- X.cs[-tt, ]
  ycv <- y.cs[-tt]
  Xtt <- X.cs[tt, ]
  ytt <- y.cs[tt]
  ids <- ids[-tt]

  # fit
  mm <- lm(ycv ~ Xcv)
  a <- mm$coeff[1]
  b <- mm$coeff[-1]

  ytt.ols <- apply(Xtt, 1, my.model, a=a,b=b,K=1)
  K.cv[j] <- abs(cov(ytt,ytt.ols))/var(ytt.ols)
  ytt.cv <- apply(Xtt, 1, my.model, a=a,b=b,K=K.cv[j])
  pmse.cv[j] <- mean((ytt - ytt.cv)^2)
}
K$crosval[i] <- K.cv[which(pmse.cv == min(pmse.cv))]


# PMSE (VS)
y.ols <- apply(X.vs, 1, my.model, a=alpha.ols,b=beta.ols,K=1)
y.copas <- apply(X.vs, 1, my.model, a=alpha.ols,b=beta.ols,K=K$copas[i])
y.crosval <- apply(X.vs, 1, my.model, a=alpha.ols,b=beta.ols,K=K$crosval[i])
PMSE$ols[i] <- mean((y.vs - y.ols)^2)
PMSE$copas[i] <- mean((y.vs - y.copas)^2)
PMSE$crosval[i] <- mean((y.vs - y.crosval)^2)
}

save(list=ls(), file='copas.Rdata')
#load('copas.Rdata')


mean(K$copas)
mean(K$crosval)

PMSE <- data.frame(PMSE)
boxplot(PMSE$copas, PMSE$ols, names=names(PMSE)[c(2,1)])
t.test(PMSE$copas, PMSE$ols, alter = 'less')
t.test(PMSE$copas, PMSE$crosval, alter = 'less')
#t.test(PMSE.copas, PMSE.crosval, alter = 'greater') #non-copas

######################
#check
######################
# randomness
#hist(X.cs)
#hist(eps)

# sig2
t.test(sig2hat, mu=sig2)

# V vs. cov
t.test(Vnorm, covXnorm)
