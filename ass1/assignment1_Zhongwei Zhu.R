setwd("/Users/aries_zzw/Documents/STAT215B/Assignment 1/")
babies <- read.table("babies.data", header=T) #Load the data as a data frame

# Prepare the data
# Figure out the missing values
attach(babies)
hist(bwt)
summary(bwt) #bwt is OK
hist(gestation) #Large value at ~1000
summary(gestation) #999 represents the missing value
babies$gestation[babies$gestation==999] <- NA
hist(parity)
summary(parity) #parity is OK
hist(age) #Large value at ~100
summary(age) #99 represents the missing value
babies$age[babies$age==99] <- NA
hist(height) #Large value at ~100
summary(height) #99 represents the missing value
babies$height[babies$height==99] <- NA
hist(weight) #Large value at ~1000
summary(weight) #999 represents the missing value
babies$weight[babies$weight==999] <- NA
hist(smoke) #Large value at ~9
summary(smoke) #9 represents the missing value
babies$smoke[babies$smoke==9] <- NA
detach(babies)
# Missing values are all replaced by NA

# Descriptive statistics or graphs
# Histograms of all variables
attach(babies)
png("Figure1.png", 1700, 800)
par(mfrow=c(2,4), mar=c(5,5,2,3), cex=1.6, cex.lab=1.2)
hist(bwt, freq=F, main="Histogram of Newborn Weights", xlab="Newborn Weight (Ounce)")
hist(gestation, freq=F, main="Histogram of Pregnancy Lengths", xlab="Pregnancy Length (Day)")
hist(as.numeric(parity), freq=F, main="Histogram of Parity", xlab="Parity")
hist(age, freq=F, main="Histogram of Mother's Ages", xlab="Mother's Age at Conception (Year)")
hist(height, freq=F, main="Histogram of Mother's Heights", xlab="Mother's Height (Inch)")
hist(weight, freq=F, main="Histogram of Mother's Weights", xlab="Mother's Weight (Pound)")
hist(as.numeric(smoke), freq=F, main="Histogram of Mother's Smoke", xlab="Mother's Smoke")
dev.off()
# Pairwise plots
panel.correlation <- function(x, y) {
	par(usr=c(0,1,0,1)) #Set the extremes of the coordinates
	a <- c(which(is.na(x)==T), which(is.na(y)==T))
	if(length(a)>0) {x <- x[-a]; y <- y[-a]}
	cor.output <- format(cor(x, y), digits=2)
	text(0.5, 0.5, cor.output, cex=4)
} #Plot the correlation values in the lower panel
panel.point <- function(x, y){
	points(x, y, pch=20, col=8, cex=2)
} #Plot the scattered plots in the upper panel
png("Figure2.png", 1200, 1200)
par(cex.axis=2)
pairs(babies, lower.panel=panel.correlation, upper.panel=panel.point, cex.labels=4, gap=1)
dev.off()

# Convert parity and smoke from numeric values to factors
babies[,3] <- as.factor(babies[,3]) #parity
babies[,7] <- as.factor(babies[,7]) #smoke
summary(babies)
detach(babies); attach(babies)

# Claim 1
# 1.1
png("Figure3.png", 1200, 600)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.6, cex.lab=1.3)
boxplot(gestation[smoke==1], gestation[smoke==0], col=c(2, 5), names=c("Smoking", "Non-Smoking"), range=0, main="Boxplot of Pregnancy Lengths", cex.axis=1.3, cex.main=1.5)
hist(gestation[smoke==1], freq=F, main="Histogram of Pregnancy Lengths", xlab="Pregnancy Length (Day)", xlim=c(140, 360), breaks=20, col="#FF000050", cex.main=1.5)
hist(gestation[smoke==0], freq=F, add=T, breaks=30, col="#00FFFF50")
legend("topleft", pch=15, col=c("#FF000080", "#00FFFF"), c("Smoking Mothers", "Non-Smoking Mothers"))
dev.off()
# 1.2
premature <- rep(0, nrow(babies)) #Create a vector for premature birth
premature[gestation<252] <- 1 #Birth prior to the 37th week (36 weeks=252 days) is premature
premature[is.na(gestation)==T] <- NA #Missing value in gestation should be denoted as NA
premature <- as.factor(premature) #Factorize the premature argument
babies <- cbind(babies, premature) #Add the factor to the data frame
table.s.p <- table(smoke, premature) #Create a talbe for comparison
dimnames(table.s.p)[[1]] <- c("No", "Yes")
dimnames(table.s.p)[[2]] <- c("No", "Yes")
# 1.3
png("Figure4.png", 800, 800)
par(cex.lab=2.4, cex.main=2.4)
mosaicplot(table.s.p, color=c(0,8), xlab="Smoking", ylab="Premature Birth", cex.axis=1.8, main="Comparison between Smoking and Premature Birth") # Mosaic plot
dev.off()
# 1.4
chisq.test(table.s.p)$p.value #Chi-squared test, p-value=0.753
fisher.test(table.s.p)$p.value #Fisher's exact test, p-value=0.683
# 1.5
gest.s <- gestation[which(smoke==1)]
gest.n <- gestation[which(smoke==0)]
mean(gest.s, na.rm=T); mean(gest.n, na.rm=T)# 278.0 for smokers and 280.2 for non-smokers 
# Check normality for t-test
png("Figure5.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,6,5,3), cex.lab=2, cex.axis=2, cex.main=2.2)
qqnorm(gest.s, pch=19, main="Q-Q Plot of Pregnancy Lengths of Smoking Mothers")
qqline(gest.s, lwd=2, col=2) #bwt of smoke mothers satify the normal distribution
qqnorm(gest.n, pch=19, main="Q-Q Plot of Pregnancy Lengths of Non-Smoking Mothers")
qqline(gest.n, lwd=2, col=2) #bwt of smoke mothers deviate from the normal distribution
dev.off()
var(gest.s, na.rm=T); var(gest.n, na.rm=T) #Variances are different for the two categories
t.test(gest.s, gest.n, alternative="less")$p.value #t-test, 0.008
wilcox.test(gest.s, gest.n, alternative="less")$p.value #wilcox-test, 0.0004
# 1.6

# Claim 2
# 2.1
bwt.s <- bwt[which(smoke==1)]
bwt.n <- bwt[which(smoke==0)]
dif.s <- mean(bwt.n) - mean(bwt.s) #Difference between avg bwt between non-smoking and smoking mothers, 8.94
bwt.f <- bwt[which(parity==1)]
bwt.nf <- bwt[which(parity==0)]
dif.f <- mean(bwt.nf) - mean(bwt.f) #Difference between avg bwt between non-first-borns and first-borns, 1.93
table.s.f <- table(smoke, parity)
var1 <- var(bwt.s)/length(bwt.s) + var(bwt.n)/length(bwt.n) + var(bwt.f)/length(bwt.f) + var(bwt.nf)/length(bwt.nf)
var2 <- 2*(-table.s.f[1,1]/(length(bwt.n)*length(bwt.nf)) + table.s.f[1,2]/(length(bwt.n)*length(bwt.f)) + table.s.f[2,1]/(length(bwt.s)*length(bwt.nf)) - table.s.f[2,2]/(length(bwt.s)*length(bwt.f)))*var(bwt)
var.s.f <- var1 + var2 #variance of differences in average difference
stat.s.f <- (dif.s-dif.f)/sqrt(var.s.f) #Wald test statistic 4.50
2*pnorm(-abs(stat.s.f)) #Two-sided p-value, 6.8E-6
# 2.2
tall <- rep(0, nrow(babies))
tall[height>median(height, na.rm=T)] <- 1
tall[is.na(height)==T] <- NA
tall <- as.factor(tall) #Divide the mothers into tall and short based on the median height
bwt.t <- bwt[which(tall==1)]
bwt.nt <- bwt[which(tall==0)]
dif.t <- mean(bwt.t) - mean(bwt.nt) #Difference between avg bwt between tall and short mothers, 5.98
table.s.t <- table(smoke, tall)
var1 <- var(bwt.s)/length(bwt.s) + var(bwt.n)/length(bwt.n) + var(bwt.t)/length(bwt.t) + var(bwt.nt)/length(bwt.nt)
var2 <- 2*(table.s.t[1,1]/(length(bwt.n)*length(bwt.nt)) - table.s.t[1,2]/(length(bwt.n)*length(bwt.t)) - table.s.t[2,1]/(length(bwt.s)*length(bwt.nt)) + table.s.t[2,2]/(length(bwt.s)*length(bwt.t)))*var(bwt)
var.s.t <- var1 + var2 #variance of differences in average difference
stat.s.t <- (dif.s-dif.t)/sqrt(var.s.t) #Wald test statistic 1.97
2*pnorm(-abs(stat.s.t)) #Two-sided p-value, 0.049
# 2.3
heavy <- rep(0, nrow(babies))
heavy[weight>median(weight, na.rm=T)] <- 1
heavy[is.na(weight)==T] <- NA
heavy <- as.factor(heavy) #Divide the mothers into heavy and light based on the median weight
bwt.h <- bwt[which(heavy==1)]
bwt.nh <- bwt[which(heavy==0)]
dif.h <- mean(bwt.h) - mean(bwt.nh) #Difference between avg bwt between heavy and light mothers, 5.28
table.s.h <- table(smoke, heavy)
var1 <- var(bwt.s)/length(bwt.s) + var(bwt.n)/length(bwt.n) + var(bwt.h)/length(bwt.h) + var(bwt.nh)/length(bwt.nh)
var2 <- 2*(table.s.h[1,1]/(length(bwt.n)*length(bwt.nh)) - table.s.h[1,2]/(length(bwt.n)*length(bwt.h)) - table.s.h[2,1]/(length(bwt.s)*length(bwt.nh)) + table.s.h[2,2]/(length(bwt.s)*length(bwt.h)))*var(bwt)
var.s.h <- var1 + var2 #variance of differences in average difference
stat.s.h <- (dif.s-dif.h)/sqrt(var.s.h) #Wald test statistic 2.53
2*pnorm(-abs(stat.s.h)) #Two-sided p-value, 0.011
# 2.4
png("Figure6.png", 1200, 1200)
par(mfrow=c(2,2), mar=c(5,5,3,3), cex=1.8, cex.lab=1.3, cex.main=1.3)
hist(bwt.s, freq=F, breaks=20, main="Histogram of Newborn Weights", xlab="Newborn Weights (Ounce)", xlim=c(50, 180), ylim=c(0, 0.035), col="#FF000050")
hist(bwt.n, freq=F, breaks=20, add=T, col="#00FFFF50")
legend("topleft", pch=15, col=c("#FF000080", "#00FFFF"), c("Smoking Mothers", "Non-Smoking Mothers"))
hist(bwt.f, freq=F, breaks=20, main="Histogram of Newborn Weights", xlab="Newborn Weights (Ounce)", xlim=c(50, 180), ylim=c(0, 0.035), col="#FF000050")
hist(bwt.nf, freq=F, breaks=20, add=T, col="#00FFFF50")
legend("topleft", pch=15, col=c("#FF000080", "#00FFFF"), c("First-Borns", "Non-First-Borns"))
hist(bwt.t, freq=F, breaks=20, main="Histogram of Newborn Weights", xlab="Newborn Weights (Ounce)", xlim=c(50, 180), ylim=c(0, 0.035), col="#FF000050")
hist(bwt.nt, freq=F, breaks=20, add=T, col="#00FFFF50")
legend("topleft", pch=15, col=c("#FF000080", "#00FFFF"), c("Tall Mothers", "Short Mothers"))
hist(bwt.h, freq=F, breaks=20, main="Histogram of Newborn Weights", xlab="Newborn Weights (Ounce)", xlim=c(50, 180), ylim=c(0, 0.035), col="#FF000050")
hist(bwt.nh, freq=F, breaks=20, add=T, col="#00FFFF50")
legend("topleft", pch=15, col=c("#FF000080", "#00FFFF"), c("Heavy Mothers", "Light Mothers"))
dev.off()
png("Figure7.png", 1200, 1200)
par(mfrow=c(2,2), mar=c(5,5,5,3), cex=1.2)
boxplot(bwt.s, bwt.n, col=c(2, 5), names=c("Smoking", "Non-Smoking"), main="Boxplot of Newborn Weights", cex.axis=1.8, cex.main=1.8, pch=19)
boxplot(bwt.f, bwt.nf, col=c(2, 5), names=c("First-Born", "Non-First-Born"), main="Boxplot of Newborn Weights", cex.axis=1.8, cex.main=1.8, pch=19)
boxplot(bwt.t, bwt.nt, col=c(2, 5), names=c("Tall", "Short"), main="Boxplot of Newborn Weights", cex.axis=1.8, cex.main=1.8, pch=19)
boxplot(bwt.h, bwt.nh, col=c(2, 5), names=c("Heavy", "Light"), main="Boxplot of Newborn Weights", cex.axis=1.8, cex.main=1.8, pch=19)
dev.off()
# 2.5
fit1 <- lm(bwt~height+weight+parity, subset=(is.na(smoke)==F))
summary(fit1)
# 2.6
fit2 <- lm(bwt~height+weight+parity+smoke)
summary(fit2)
# Plot the residuals over fitted values and make the Q-Q plots for residuals
png("Figure8.png", 1500, 800)
par(mfrow=c(1,2), mar=c(5,5,3,3), cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.8)
plot(x=fit1$fitted.values, y=fit1$residuals, pch=20, xlab="Fitted Values", ylab="Residuals", main="Residuals of Fit without Smoking")
abline(h=2*c(sd(fit1$residuals), -sd(fit1$residuals)), lwd=2, lty=2)
plot(x=fit2$fitted.values, y=fit2$residuals, pch=20, xlab="Fitted Values", ylab="Residuals", main="Residuals of Fit with Smoking")
abline(h=2*c(sd(fit2$residuals), -sd(fit2$residuals)), lwd=2, lty=2)
dev.off()
# ANOVA
anova(fit1)
anova(fit2)
anova(fit1, fit2)