setwd("/Users/aries_zzw/Documents/STAT215B/Final")
rm(list=ls())

# Read the Dataset
ibm <- read.table("IBMC.csv", sep=",", header=T)
ibm <- ibm[,-1]
id <- ibm$lot_ID*100 + ibm$wafer_ID # Define a new id for each wafer
wafer.id <- unique(id)
die.loc <- unique(ibm$chipX*100 + ibm$chipY) # Define new positions for each die
wafer.amount <- length(wafer.id) # Total amount of wafers

# Ture fault wafers, expert label
outlierTrue <- c(421,507,523,707,1605,1607,2201,2205,508,703,421,1305,1412,1507,1604,1606,1207,1201,607,514,511,509,505,1909,1915,1705,517,520,521,1217,1311,1415,1901,1906,1911,503,1312,1413,515)
outlierTrue <- sort(unique(outlierTrue))

# Use die mean for the analysis
ibm.data <- as.matrix(cbind(id, ibm[,-c(1:2)]))
die.avg <- rowMeans(ibm.data[,4:ncol(ibm.data)])
die.sd <- apply(ibm.data[,4:ncol(ibm.data)], 1, sd)
plot(die.avg, die.sd, pch=19)
wafer.avg <- ave(die.avg, as.factor(ibm.data[,1]))
w.d.res <- matrix(die.avg-wafer.avg, wafer.amount, byrow=T, dimnames=list(wafer.id, die.loc))
# Each element in w.d.res represents the mean of 14 measures on a given die. Each row and column represents a wafer and a certain position.

# PCA on die level
die.pca <- prcomp(ibm.data[,-(1:3)], T, T, F)
cumsum(die.pca$sdev^2)/sum(die.pca$sdev^2) #First component takes over 85% of the variance
die.pca$rot[,1] #Coefficients are similar for each measurement, so using mean is valid

# Divide data into training and test sets
# train.size <- round(wafer.amount*0.8)
# test.size <- wafer.amount-train.size
# train.id <- sample(wafer.id, train.size)
# train.id <- train.id[order(train.id)]
# test.id <- wafer.id[wafer.id %in% train.id == F]
# ibm.train <- ibm.data[id %in% train.id,]
# ibm.test <- ibm.data[id %in% test.id,]

# PCA on the die residuals
pca.r <- prcomp(w.d.res, T, T ,F)
cumsum(pca.r$sdev^2)/sum(pca.r$sdev^2) # First 3 components contain 62% of the total variance
label.pca.r <- rep(F, wafer.amount)
require(sfsmisc)
get.fault <- function(i, x, data, bw){
	data <- data[,i]
	hist.data <- hist(data, breaks=20, plot=F)
	dens.data <- density(data, bw=bw)
	m <- seq(round(min(data), 2), round(max(data), 2), 0.01)
	area <- sapply(m, function(a) integrate.xy(dens.data$x, dens.data$y, b=a))
	if(max(area)<(1-x/2)) ub<- round(max(data), 2) else ub <- m[which.max(area>=(1-x/2))]
	lb <- m[which.min(area<=(x/2))]
	return(c("lower"=lb, "upper"=ub))
}
bound <- sapply(1:3, get.fault, x=0.05, data=pca.r$x, bw=0.2)
label.pca.r[pca.r$x[,1]>bound[2,1] | pca.r$x[,1]<bound[1,1]] <- T
label.pca.r[pca.r$x[,2]>bound[2,2] | pca.r$x[,2]<bound[1,2]] <- T
label.pca.r[pca.r$x[,3]>bound[2,3] | pca.r$x[,3]<bound[1,3]] <- T
png("Figure1.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(pca.r$x[,1], pca.r$x[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component", col=label.pca.r+1)
abline(h=bound[,2], v=bound[,1], lwd=2, lty=2)
plot(pca.r$x[,2], pca.r$x[,3], pch=19, xlab="Second Principal Component", ylab="Third Principal Component", col=label.pca.r+1)
abline(h=bound[,3], v=bound[,2], lwd=2, lty=2)
dev.off()
get.sum <- function(x){
  tot <- sum(x) # Number of fault wafers identified
  FR <- sum((wafer.id[x==T] %in% outlierTrue) == F) # False rate
  mis <- sum((outlierTrue %in% wafer.id[x==T]) == F) # Misclassification rate
  return(list(c("tot"=tot, "FR"=FR, "mis"=mis), wafer.id[x==T]))
}
get.sum(label.pca.r)
num <- function(x, data, bw, fn){
	bound <- sapply(1:3, fn, x=x, data=data, bw=bw)
	label <- rep(F, wafer.amount)
	label[data[,1]>bound[2,1] | data[,1]<bound[1,1]] <- T
	label[data[,2]>bound[2,2] | data[,2]<bound[1,2]] <- T
	label[data[,3]>bound[2,3] | data[,3]<bound[1,3]] <- T
	result <- get.sum(label)
	return(result[[1]])
}
fdr.r <- sapply(0.01*(0:100), num, data=pca.r$x, bw=0.2, fn=get.fault)
bound.1 <- get.fault(1, 0.05, pca.r$x, 0.1)
png("Figure2.png", 1000, 800)
par(mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
hist(pca.r$x[,1], breaks=40, freq=F, xlab="First Principal Component", main="Histogram of the First Component")
lines(density(pca.r$x[,1], bw=0.1), col=2, lwd=2)
abline(v=bound.1, col=4, lwd=2)
dev.off()
# Estimate variograms
require(geoR)
geo.data <- cbind(ibm.data[,1:3], die.avg)
geo.fit <- sapply(1:wafer.amount, function(i){
  gR <- geo.data[geo.data[,1] == wafer.id[i],]
  gR <- as.geodata(gR, coords.col=2:3, data.col=4)
  vario <- variog(gR, estimate="modulus")
  return(vario)
})
geo.fit <- t(geo.fit); rownames(geo.fit) <- wafer.id
geo.u <- sapply(1:wafer.amount, function(i) geo.fit[i,1][[1]])
sum(geo.u!=rowMeans(geo.u)) # All u's are the same for all the wafers
geo.v <- sapply(1:wafer.amount, function(i) geo.fit[i,2][[1]])
geo.v <- t(geo.v); rownames(geo.v) <- wafer.id

# PCA on variograms
pca.v <- prcomp(geo.v, T, T, F)
cumsum(pca.v$sdev^2)/sum(pca.v$sdev^2) # First 3 components contain 98% of the total variance
label.pca.v <- rep(F, wafer.amount)
bound <- sapply(1:3, get.fault, x=0.05, data=pca.v$x, bw=0.05)
label.pca.v[pca.v$x[,1]>bound[2,1] | pca.v$x[,1]<bound[1,1]] <- T
label.pca.v[pca.v$x[,2]>bound[2,2] | pca.v$x[,2]<bound[1,2]] <- T
label.pca.v[pca.v$x[,3]>bound[2,3] | pca.v$x[,3]<bound[1,3]] <- T
png("Figure3.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(pca.v$x[,1], pca.v$x[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component", col=label.pca.v+1)
abline(h=bound[,2], v=bound[,1], lwd=2, lty=2)
plot(pca.v$x[,2], pca.v$x[,3], pch=19, xlab="Second Principal Component", ylab="Third Principal Component", col=label.pca.v+1)
abline(h=bound[,3], v=bound[,2], lwd=2, lty=2)
dev.off()
get.sum(label.pca.v)
fdr.v <- sapply(0.01*(0:100), num, data=pca.v$x, bw=0.05, fn=get.fault)

# Kernel PCA
require(kernlab)
# Die residuals, normal kernel
kpca.r <- kpca(w.d.res, kpar=list(sigma=0.5))
rot.k.r <- rotated(kpca.r)
k.kpca.r <- kmeans(rot.k.r, 2)
rot.k.r2 <- rotated(kpca(w.d.res, kpar=list(sigma=1)))
k.kpca.r2 <- kmeans(rot.k.r2, 2)
label.kpca.r <- label.kpca.r2 <- rep(F, wafer.amount)
label.kpca.r[k.kpca.r$cluster==1] <- T
label.kpca.r2[k.kpca.r2$cluster==1] <- T
png("Figure4.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(rot.k.r[,1], rot.k.r[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component", col=label.kpca.r+1)
plot(rot.k.r2[,1], rot.k.r2[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component", col=label.kpca.r2+1)
dev.off()
get.sum(label.kpca.r); get.sum(label.kpca.r2)

# Variogram, normal kernel
kpca.v <- kpca(geo.v, kpar=list(sigma=1))
rot.k.v <- rotated(kpca.v)
label.kpca.v <- rep(F, wafer.amount)
get.fault.k <- function(i, x, data, bw){
	data <- data[,i]
	hist.data <- hist(data, breaks=20, plot=F)
	dens.data <- density(data, bw=bw)
	m <- seq(round(min(data), 2), round(max(data), 2), 0.01)
	area <- sapply(m, function(a) integrate.xy(dens.data$x, dens.data$y, b=a))
	if(max(area)<(1-x)) ub<- round(max(data), 2) else ub <- m[which.max(area>=(1-x))]
	return(ub)
}
bound <- sapply(1:2, get.fault.k, x=0.05, data=rot.k.v, bw=0.5)
label.kpca.v[rot.k.v[,1]>bound[1]] <- T
label.kpca.v[rot.k.v[,2]>bound[2]] <- T
png("Figure5.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(rot.k.v[,1], rot.k.v[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component",  col=label.kpca.v+1)
abline(h=bound[2], v=bound[1], lwd=2, lty=2)
plot(rot.k.v[,2], rot.k.v[,3], pch=19, xlab="Second Principal Component", ylab="Third Principal Component",  col=label.kpca.v+1)
abline(v=bound[2], lwd=2, lty=2)
dev.off()
get.sum(label.kpca.v)
num.k <- function(x, data, bw, fn){
	bound <- sapply(1:2, fn, x=x, data=data, bw=bw)
	label <- rep(F, wafer.amount)
	label[data[,1]>bound[1]] <- T
	label[data[,2]>bound[2]] <- T
	result <- get.sum(label)
	return(result[[1]])
}
fdr.k.v <- sapply(0.01*(0:100), num.k, data=rot.k.v, bw=0.5, fn=get.fault.k)
png("Figure6.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(fdr.r[3,], fdr.r[2,], type="s", xlab="Number of Misclassifications", ylab="Number of False Alarms", lwd=3)
lines(fdr.v[3,], fdr.v[2,], type="s", lwd=3, col=2)
legend("topright", c("Die.residual", "Variogram"), col=c(1,2), lty=1, lwd=2, cex=1.5, bty="n")
plot(fdr.v[3,], fdr.v[2,], type="s", xlab="Number of Misclassifications", ylab="Number of False Alarms", lwd=3, col=2)
lines(fdr.k.v[3,], fdr.k.v[2,], type="s", lwd=3, col=4)
legend("topright", c("Linear PCA", "Kernel PCA"), col=c(2,4), lty=1, lwd=2, cex=1.5, bty="n")
dev.off()

label.true <- rep(0, wafer.amount)
label.true[wafer.id%in%outlierTrue] <- 1
