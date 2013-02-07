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
head(babies)
#hist(babies$gestation)

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
head(babies)



#### Claim 1: gestation ~ smoke
## 1. graphical comparison: smoke vs. non-smoke
attach(babies)
typ="h";lwd=6;col="blue"
xlab="Days of pregnancy"; title="Gestation Distribution"
#par(mfrow=c(2,1),cex.main= 1)
#hist(gestation[smoke == "No"], type=typ,col=col,lwd=lwd, main=title,xlab=xlab,ylab="Non-smoking mother")
#hist(gestation[smoke == "Yes"], type=typ,col=col,lwd=lwd, main=title,xlab=xlab,ylab="Smoking mother")


## 2. tabular comparison
tb.PS <- table(premature, smoke)
tb.PS


## 3. visualization of 2.


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
##??
t.test(bwt ~ smoke)
t.test(bwt ~ parity)


## 2. bwt difference ~ smoke/height
t.test(bwt ~ smoke)
t.test(bwt ~ f.height)



## 3. bwt difference ~ smoke/weight
t.test(bwt ~ smoke)
t.test(bwt ~ f.weight)



## 4. visual comparison


## 5. LS (without smoking)
fit1 = lm(bwt ~ height + weight + parity, na.action= na.exclude)
summary(fit1)



## 6. LS (with smoking)
fit2 = lm(bwt ~ height + weight + parity + smoke, na.action= na.exclude)
summary(fit2)


# Residual-plot, important tool in checking assumptions of model.
plot(fitted(fit2), resid(fit2))
# Normality-assumption of residuals reasonable?
hist(resid(fit2))
shapiro.test(fit2$resid)

anova(fit1)
anova(fit2)

## 7. multivariate vs. univariate


## 8. other statistics



## 9. ggplot2


detach(babies)









