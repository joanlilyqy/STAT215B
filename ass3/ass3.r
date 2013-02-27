#### Ying Qiao
#### SID: 21412301
#### STAT 215B, Assignment 3, 3/1/2013

rm(list=ls())

#### Detective Work
## N : total subject counts for various subgroups/margins
unsub <- read.table("part6_907.txt", header=T)














# 1. random generation for the Weibull distribution with
# parameters shape (b) and scale (a) using the inverse CDF 
# G^{-1}_{a,b}(x) = a{-log(1-x)}^{1/b}, a>0 b>0
# Invalid arguments will result in return value 'NaN', with a warning.

my.rweibull <- function (n, shape, scale = 1){
  if (n <= 0 || shape <= 0 || scale <= 0){
    warning("Invalid arguments! All arguments should be positive.")
    return ("NaN")
  }

  u <- runif(n)
  x <- scale * (-log(1 - u))^(1/shape)
  return (x)
}


# 2. creates a Kaplan-Meier survival-function estimate
# S(t_i) = {(N(t_i)-E(t_i))/N(t_i)}*S(t_{i-1})
# input: event times; cencoring indicator (TRUE for death, FALSE for censor)
# output: estimate points t, survival time estimates s, estimates survival function my.surv(eval)

my.KMsurv <- function (times, censor){
  if (length(times) <= 0 || length(censor) <= 0)
    stop("Invalid input vectors!")
  if (length(times) != length(censor))
    stop("Mismatch of dimension in input vectors!")

  t <- c(0, sort(unique(times[censor]))) # prepare t estimates points
  s <- rep(0, length(t)) # create empty survival time estimates vector
  # K-M method
  s[1] <- 1 # initialize 
  for (i in 2:length(t)){
    n <- length(times[times >= t[i]]) # num of subjects at risk at time t
    e <- length(times[times == t[i]]) # num of events(deaths) at time t
    s[i] <- (n-e)/n * s[i-1]
  }
  
  # estimated survival function S(t)
  my.surv <- function (eval){
    return (s[t == eval])
  }

  # return a list of related KM objects
  km <- list(t = t, surv = s, s.func = my.surv)
  return (km)
}


## Simulation study setup
# N=1000 | 1/2
# X ~ G_{a=3,b=2} -- treatment
# Y ~ G_{a=2,b=2} -- control
set.seed(0)
X <- my.rweibull(500,shape=2,scale=3)
set.seed(100)
Y <- rweibull(500,shape=2,scale=2)
cen <- rep(TRUE, 500)
# true curve
x.km0 <- my.KMsurv(X, cen)
y.km0 <- my.KMsurv(Y, cen)
### R built-in function check
### require(survival)
### kmfit <- survfit(Surv(X, cen) ~ 1)



# 3. simple simulation censoring: cut-off at 5 years
cen.x <- X <= 5
cen.y <- Y <= 5
x.km1 <- my.KMsurv(X, cen.x)
y.km1 <- my.KMsurv(Y, cen.y)
# visual comparisons
par(mfrow=c(1,2),cex.main=1)
typ='s'; lwd=2; col0='blue';lty0=1;col='red';lty=4
xlab="Time (year)"; ylab="Survival Rate"; ylim=c(0,1)
plot(x.km0$t, x.km0$surv, type=typ,col=col0,lwd=lwd,lty=lty0, xlim=c(0,max(X)),ylim=ylim, xlab=xlab,ylab=ylab, main="Survival Curves for Treatment Group X")
lines(x.km1$t, x.km1$surv, type=typ,col=col,lwd=lwd,lty=lty)
legend("topright",c("True","Estimated"),col=c(col0,col),lty=c(lty0,lty))
plot(y.km0$t, y.km0$surv, type=typ,col=col0,lwd=lwd,lty=lty0, xlim=c(0,max(X)),ylim=ylim, xlab=xlab,ylab=ylab, main="Survival Curves for Control Group Y")
lines(y.km1$t, y.km1$surv, type=typ,col=col,lwd=lwd,lty=lty)
legend("topright",c("True","Estimated"),col=c(col0,col),lty=c(lty0,lty))


# 4. advanced simulation censoring: cut-off at iid Zi~Exp(mean=10)
set.seed(0);   Z <- rexp(500, rate=1/10)
cen.x <- X <= Z
set.seed(100); Z <- rexp(500, rate=1/10)
cen.y <- Y <= Z
x.km2 <- my.KMsurv(X, cen.x)
y.km2 <- my.KMsurv(Y, cen.y)
# visual comparisons
par(mfrow=c(1,2),cex.main=1)
plot(x.km0$t, x.km0$surv, type=typ,col=col0,lwd=lwd,lty=lty0, xlim=c(0,max(X)),ylim=ylim, xlab=xlab,ylab=ylab, main="Survival Curves for Treatment Group X")
lines(x.km2$t, x.km2$surv, type=typ,col=col,lwd=lwd,lty=lty)
legend("topright",c("True","Estimated with Z"),col=c(col0,col),lty=c(lty0,lty))
plot(y.km0$t, y.km0$surv, type=typ,col=col0,lwd=lwd,lty=lty0, xlim=c(0,max(X)),ylim=ylim, xlab=xlab,ylab=ylab, main="Survival Curves for Control Group Y")
lines(y.km2$t, y.km2$surv, type=typ,col=col,lwd=lwd,lty=lty)
legend("topright",c("True","Estimated with Z"),col=c(col0,col),lty=c(lty0,lty))


# 5. advanced simulation censoring: cut-off at iid Zi|Xi<2~Exp(mean=10); Zi|Xi>=2~Exp(mean=5)
Z <- rep(0, 500)
cut.x <- X < 2; n.l2x <- length(X[cut.x])
set.seed(0);   Z[cut.x]  <- rexp(n.l2x, rate=1/10)
set.seed(10);  Z[!cut.x] <- rexp(500-n.l2x, rate=1/5)
cen.x <- X <= Z
Z <- rep(0, 500)
cut.y <- Y < 2; n.l2y <- length(Y[cut.y])
set.seed(100); Z[cut.y]  <- rexp(n.l2y, rate=1/10) 
set.seed(110); Z[!cut.y] <- rexp(500-n.l2y, rate=1/5)
cen.y <- Y <= Z
x.km3 <- my.KMsurv(X, cen.x)
y.km3 <- my.KMsurv(Y, cen.y)
# visual comparisons
par(mfrow=c(1,2),cex.main=1)
plot(x.km0$t, x.km0$surv, type=typ,col=col0,lwd=lwd,lty=lty0, xlim=c(0,max(X)),ylim=ylim, xlab=xlab,ylab=ylab, main="Survival Curves for Treatment Group X")
lines(x.km3$t, x.km3$surv, type=typ,col=col,lwd=lwd,lty=lty)
legend("topright",c("True","Estimated with Z|X"),col=c(col0,col),lty=c(lty0,lty))
plot(y.km0$t, y.km0$surv, type=typ,col=col0,lwd=lwd,lty=lty0, xlim=c(0,max(X)),ylim=ylim, xlab=xlab,ylab=ylab, main="Survival Curves for Control Group Y")
lines(y.km3$t, y.km3$surv, type=typ,col=col,lwd=lwd,lty=lty)
legend("topright",c("True","Estimated with Z|Y"),col=c(col0,col),lty=c(lty0,lty))
