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
A.hat <- sum(mu^2) / (n-2)
TSE <- (n*A.hat+2)/(A.hat+1) #8.05 vs. 8.13

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

#####################################################################
rm(list=ls())

#### Shrinking radon
## load data
srrs <- read.table("srrs2.dat", sep=",", header=T)
#head(srrs)

## cleaning
base <- srrs[srrs$floor == 0, ]
mn <- base[base$state == 'MN', ]
obs.ct <- table(mn$county)
maj.ct <- obs.ct[obs.ct >= 10]
N <- length(maj.ct)
radon <- mn[mn$county %in% names(maj.ct), ]
dim(radon)
#head(radon)

## split set
n <- 5
train <- c()
for (i in 1:N)
{
tmp <- radon[radon$county == names(maj.ct)[i], ]
set.seed(i*10)
ids <- sample(1:nrow(tmp), n)
train <- rbind(train, tmp[ids, ])
}
dim(train)
test <- radon[!(radon$idnum %in% train$idnum), ]
dim(test)

## "true" parameter value
mu <- rep(0, N)
names(mu) <- names(maj.ct)
for (i in 1:N)
{
tmp <- test[test$county == names(maj.ct)[i], ]
if ( nrow(tmp) == maj.ct[i] - n )
  mu[i] <- mean(tmp$activity)
}





## MLE estimator
mu.MLE <- rep(0, N)
wc.se <- rep(0, N) # within-county sum of squared residuals
names(mu.MLE) <- names(maj.ct)
for (i in 1:N)
{
tmp <- train[train$county == names(maj.ct)[i], ]
if ( nrow(tmp) == n ){
  mu.MLE[i] <- mean(tmp$activity)
  wc.se[i] <- sum((tmp$activity - mu.MLE[i])^2)
}
}
mu.MLE
TSE.MLE <- sum((mu.MLE - mu)^2)
TSE.MLE


## brief note
# z_{i,j} ~ N(u_i, t^2) i=1...17, j=1...5 all independent
# uMLE = z_{i,.}=avg_j{z_{i,j}}
# z_{i,.} ~ N(u_i, t^2/n_j)

## JS estimator
mu.bar <- mean(mu.MLE)
tau2 <- sum(wc.se) / (N*(n-1))
SE <- tau2 / n #note
S <- sum((mu.MLE - mu.bar)^2)
shr <- 1 - (N - 3)*SE/S

mu.JS <- rep(0, N)
names(mu.JS) <- names(maj.ct)
for (i in 1:N)
{
tmp <- train[train$county == names(maj.ct)[i], ]
if ( nrow(tmp) == n ){
  mu.JS[i] = mu.bar + shr*(mu.MLE[i] - mu.bar)
}
}
mu.JS
TSE.JS <- sum((mu.JS - mu)^2)
TSE.JS


