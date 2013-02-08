#### Ying Qiao
#### SID: 21412301
#### STAT 215B, Assignment 1, 2/8/2013

rm(list=ls())

# load and pre-process data
babies <- read.table("babies.data", header=T)
head(babies)
# replace NA
babies$gestation[babies$gestation == 999] <- as.numeric("NA")
babies$age[babies$age == 99] <- as.numeric("NA")
babies$height[babies$height == 99] <- as.numeric("NA")
babies$weight[babies$weight == 999] <- as.numeric("NA")
babies$smoke[babies$smoke == 9] <- as.numeric("NA")
# categorical var
babies$parity <- factor(babies$parity)
levels(babies$parity) <- c("not.FB","is.FB") # first-born
babies$smoke <- factor(babies$smoke)
levels(babies$smoke) <- c("No","Yes") # smokes
# check initial statistics
head(babies, n = 10)

## c1.2
babies$premature <- factor(babies$gestation < 37*7)
levels(babies$premature) <- c("not.PM","is.PM") # premature factor

## c2.2
babies$f.height <- factor(babies$height > median(babies$height, na.rm=T))
levels(babies$f.height) <- c("short","tall") # height factor

## c2.3
babies$f.weight <- factor(babies$weight > median(babies$weight, na.rm=T))
levels(babies$f.weight) <- c("light","heavy") # height factor
## check
head(babies, n = 10)
hist(babies$gestation,col="blue",lwd=2, main="Distribution of overall gestation",xlab="Days of pregnancy",ylab="Counts of mothers")


#### Claim 1: gestation ~ smoke
## 1. graphical comparison: smoke vs. non-smoke
attach(babies)
xlim=range(gestation, na.rm=T); ylim=c(0,0.04);
lwd=2;col="blue"; xlab="Days of pregnancy"; breaks=seq(140,360,10)
#boxplot
boxplot(gestation~smoke, data=babies, col=col,lwd=lwd, xlab="Smoke", ylab=xlab)
#histogram
par(mfrow=c(2,1),cex.main= 1)
hist(gestation[smoke == "No"],breaks=breaks,freq=F, xlim=xlim,ylim=ylim, col=col,lwd=lwd,main="",xlab=xlab,ylab="Non-smoking mother")
hist(gestation[smoke == "Yes"],breaks=breaks,freq=F, xlim=xlim,ylim=ylim, col=col,lwd=lwd,main="",xlab=xlab,ylab="Smoking mother")

## 2. tabular comparison
tb.PS <- table(premature, smoke)
ftable(tb.PS)

## 3. visualization of 2.
require(stats)
mosaicplot(smoke ~ premature, data=babies, main = "Gestation distribution", color = TRUE)

## 4. hypothesis test:
## H0: smoking and non-smoking mothers have the same rate of premature delivery
chisq.test(tb.PS)
fisher.test(tb.PS)

## 5. hypothesis test:
## H0: smoking and non-smoking mothers have the same overall average gestation time
## HA: smoking mothers have shorter overall average gestation time
t.test(gestation ~ smoke, alternative = "less")
wilcox.test(gestation ~ smoke, alternative = "less")

## 6. other support for claim 1





#### Claim 2: bwt ~ smoke (stronger covariate)
## 1. bwt difference ~ smoke/first-born
## assume bwt i.i.d, initial obervation
t.test(bwt ~ smoke)
t.test(bwt ~ parity)
t.test(bwt ~ f.height)
t.test(bwt ~ f.weight)

## H0: the difference in the average bwt between smoking/non-smoking mothers is the same as that between firt-borns/non-first-borns
dA = mean(bwt[smoke=="No"], na.rm=T) - mean(bwt[smoke=="Yes"], na.rm=T)
dB = mean(bwt[parity=="not.FB"], na.rm=T) - mean(bwt[parity=="is.FB"], na.rm=T)
var.dA = var(bwt[smoke=="No"], na.rm=T) + var(bwt[smoke=="Yes"], na.rm=T)
var.dB = var(bwt[parity=="not.FB"], na.rm=T) + var(bwt[parity=="is.FB"], na.rm=T)
cov.dAdB = sum(var(bwt[smoke=="No" && parity=="not.FB"], na.rm=T), var(bwt[smoke=="No" && parity=="is.FB"], na.rm=T),
               var(bwt[smoke=="Yes" && parity=="not.FB"], na.rm=T), var(bwt[smoke=="Yes" && parity=="is.FB"], na.rm=T), na.rm=T)
se.dAdB = sqrt(var.dA + var.dB - 2*cov.dAdB)

wald.test1 <- (abs(dA)-abs(dB))/se.dAdB
p.1 <- pnorm(wald.test1, lower.tail=F)

## 2. bwt difference ~ smoke/height
## H0: the difference in the average bwt between smoking/non-smoking mothers is the same as that between tall/short
dA = mean(bwt[smoke=="No"], na.rm=T) - mean(bwt[smoke=="Yes"], na.rm=T)
dB = mean(bwt[f.height=="short"], na.rm=T) - mean(bwt[f.height=="tall"], na.rm=T)
var.dA = var(bwt[smoke=="No"], na.rm=T) + var(bwt[smoke=="Yes"], na.rm=T)
var.dB = var(bwt[f.height=="short"], na.rm=T) + var(bwt[f.height=="tall"], na.rm=T)
cov.dAdB = sum(var(bwt[smoke=="No" && f.height=="short"], na.rm=T), var(bwt[smoke=="No" && f.height=="tall"], na.rm=T),
               var(bwt[smoke=="Yes" && f.height=="short"], na.rm=T), var(bwt[smoke=="Yes" && f.height=="tall"], na.rm=T), na.rm=T)
se.dAdB = sqrt(var.dA + var.dB - 2*cov.dAdB)

wald.test2 <- (abs(dA)-abs(dB))/se.dAdB
p.2 <- pnorm(wald.test2, lower.tail=F) 

## 3. bwt difference ~ smoke/weight
## H0: the difference in the average bwt between smoking/non-smoking mothers is the same as that between light/heavy
dA = mean(bwt[smoke=="No"], na.rm=T) - mean(bwt[smoke=="Yes"], na.rm=T)
dB = mean(bwt[f.weight=="light"], na.rm=T) - mean(bwt[f.weight=="heavy"], na.rm=T)
var.dA = var(bwt[smoke=="No"], na.rm=T) + var(bwt[smoke=="Yes"], na.rm=T)
var.dB = var(bwt[f.weight=="light"], na.rm=T) + var(bwt[f.weight=="heavy"], na.rm=T)
cov.dAdB = sum(var(bwt[smoke=="No" && f.weight=="light"], na.rm=T), var(bwt[smoke=="No" && f.weight=="heavy"], na.rm=T),
               var(bwt[smoke=="Yes" && f.weight=="light"], na.rm=T), var(bwt[smoke=="Yes" && f.weight=="heavy"], na.rm=T), na.rm=T)
se.dAdB = sqrt(var.dA + var.dB - 2*cov.dAdB)

wald.test3 <- (abs(dA)-abs(dB))/se.dAdB
p.3 <- pnorm(wald.test3, lower.tail=F) 


## 4. visual comparisons
par(mfrow=c(2,2),cex.main= 1)
ylab="Birth weight";ylim=range(bwt,na.rm=T)
boxplot(bwt~smoke, data=babies, ylim=ylim, col=col,lwd=lwd, xlab="Smoke", ylab=ylab)
boxplot(bwt~parity, data=babies, ylim=ylim, col=col,lwd=lwd, xlab="First-born", ylab=ylab)
boxplot(bwt~f.height, data=babies, ylim=ylim, col=col,lwd=lwd, xlab="Height", ylab=ylab)
boxplot(bwt~f.weight, data=babies, ylim=ylim, col=col,lwd=lwd, xlab="Weight", ylab=ylab)


## 5. LS (without smoking)
fit1 = lm(bwt ~ height + weight + parity, na.action= na.exclude)
# check the fit
par(mfrow=c(1,2),cex.main= 1)
# Residual-plot, important tool in checking assumptions of model.
plot(fitted(fit1), resid(fit1), type="p", col=col, lwd=lwd, xlab="Fitted values of Model 1", ylab="Residuals of Model 1", main="Fitting check")
# Normality-assumption of residuals reasonable?
lim=max(resid(fit1), na.rm=T)
hist(resid(fit1),breaks=20,freq=F,xlim=c(-lim,lim),ylim=c(0,0.025), col=col,lwd=lwd,main="Normality check",xlab="Residual values of Model 1",ylab="Density")


## 6. LS (with smoking)
fit2 = lm(bwt ~ height + weight + parity + smoke, na.action= na.exclude)

summary(fit1)
summary(fit2)

anova(fit1)
anova(fit2)

## 7. multivariate vs. univariate


## 8. other statistics


## (extra) 9. ggplot2


detach(babies)









