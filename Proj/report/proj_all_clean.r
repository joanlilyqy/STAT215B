#### STAT215B Project
#### Haotian Liu, Ying Qiao, Zhongwei Zhu
require(ggplot2)
require(geoR)
require(e1071)
require(kernlab)


######## Missing Value Interpolation ########
rm(list=ls())
ibm <- read.table("IBMBigData.csv", sep=",", header=T) #load original data
id <- ibm$lot_ID*100 + ibm$wafer_ID # Define a new id for each wafer
wid <- unique(id) # unique wafer id list
meas <- grep("psro", names(ibm)) # measurement indices
lw.id <- function(ii){# aux function : recover lot_ID and wafer_ID
  lid <- floor(ii/100)
  wid <- ii - lid*100
  return (list(lid=lid,wid=wid))
}
jpeg(file='wafer_counts.jpg')
plot(table(floor(unique(id)/100)),col='blue',lwd=2,typ='h',xlab='Lot ID',ylab='Num. of wafers within lot')
dev.off()
#### plot wafermap with missing values
pid <- lw.id(107)
wwData <- ibm[ibm$lot == pid$l & ibm$wafer == pid$w,  ]
dieMean<-rowMeans(wwData[  ,grep("psro", names(ibm))])
dm <- cbind(wwData[  ,c('chipX','chipY')], dieMean)
plotWmap <- function(dm, id){# aux function : wafer meas heatmap
  title <- sprintf("Lot%dWafer%d", id$l, id$w)
  fn <- sprintf("%s.jpeg",title)
  fig <- ggplot(dm, aes(x = chipX, y = chipY)) + scale_x_continuous(breaks=1:13) + scale_y_continuous(breaks=1:13)
  fig <- fig + geom_point(aes(colour = dieMean),shape = 15,size = 10) + ggtitle(title)
  fig <- fig + guides(colour = guide_colorbar(barwidth = 0.5, barheight = 20))
  ggsave(filename=fn, plot=fig)
}
plotWmap(dm, pid)
#### get the full wafer die measurement location
dieTT <- table(id)
notNA <- as.numeric(names(dieTT)[which(dieTT == max(dieTT))])
nna <- lw.id(notNA[2])
nna.x <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipX']
nna.y <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipY']
#### prepare full wafer dataset
wfgr <- data.frame(cbind(nna.x, nna.y))
names(wfgr) <- c('chipX','chipY')
ndie <- length(nna.x)
IBMC <- data.frame(matrix(0, nrow=length(wid)*ndie, ncol=ncol(ibm)))
names(IBMC) <- names(ibm)
#### Variogram model fit from empirical values
rng <- min(max(nna.x),max(nna.y)) # range for variogram model
nug <- 0 # nugget for variogram model
sill <- matrix(0, nrow=length(unique(id)), ncol=length(meas)) # estimated sill for variogram model
for (i in 1:length(wid)){
  for (j in 1:length(meas)){
    gid <- lw.id(wid[i])
    gR <- ibm[ibm$lot == gid$l & ibm$wafer == gid$w,  ]
    gR <- as.geodata(gR, coords.col=c('chipX','chipY'), data.col=meas[j])
    vario <- variog(gR, estimate='modulus')
    var.mod <- variofit(vario, c(max(vario$v),max(vario$u)), cov.model='cubic')
    sill[i, j] <- var.mod$cov.pars[1]
  }
}
gid <- lw.id(107)
fn <- sprintf('Lot%dWafer%d.jpeg',gid$l,gid$w)
jpeg(file=fn)
gR <- ibm[ibm$lot == gid$l & ibm$wafer == gid$w,  ]
gR <- as.geodata(gR, coords.col=c('chipX','chipY'), data.col=meas[1])
vario <- variog(gR, estimate='modulus')
var.mod <- variofit(vario, c(max(vario$v),max(vario$u)), cov.model='cubic')
plot(vario,main='Variogram');lines(var.mod)
dev.off()
#### Ordinary kriging for interpolation
for (i in 1:length(wid)){
  gid <- lw.id(wid[i])
  wR <- ibm[ibm$lot == gid$l & ibm$wafer == gid$w,  ]
  rowl <- (1:ndie) + (i-1)*ndie
  IBMC[rowl, 'lot_ID'] <- gid$l
  IBMC[rowl, 'wafer_ID'] <- gid$w
  IBMC[rowl, c('chipX','chipY')] <- wfgr
  for (j in 1:length(meas)){    
    gR <- as.geodata(wR, coords.col=c('chipX','chipY'), data.col=meas[j])
    ctrl <- krige.control(type.krige='ok',cov.model='cubic',
                          cov.pars=c(sill[i,j],rng),nugget=nug)
    krg.mod <- krige.conv(gR, locations=wfgr, krige=ctrl)
    intpol <- krg.mod$pred
    IBMC[rowl, meas[j]] <- intpol
  }
}
write.csv(IBMC, file='IBMC.csv')
#### plot wafermap with complete values
pid <- lw.id(107)
wwData <- IBMC[IBMC$lot == pid$l & IBMC$wafer == pid$w,  ]
dieMean<-rowMeans(wwData[  ,grep("psro", names(IBMC))])
dm <- cbind(wwData[  ,c('chipX','chipY')], dieMean)
plotWmap(dm, pid)


######## PCA and Kernel PCA ########
rm(list=ls())
#### Reload the "complete" dataset
ibm <- read.table("IBMC.csv", sep=",", header=T)
ibm <- ibm[,-1]
id <- ibm$lot_ID*100 + ibm$wafer_ID # Define a new id for each wafer
wafer.id <- unique(id) # Wafer id list
die.loc <- unique(ibm$chipX*100 + ibm$chipY) # Define new positions for each die
wafer.amount <- length(wafer.id) # Total amount of wafers
#### Ture fault wafers, expert label
outlierTrue <- c(421,507,523,707,1605,1607,2201,2205,508,703,421,1305,1412,1507,
                 1604,1606,1207,1201,607,514,511,509,505,1909,1915,1705,517,520,
                 521,1217,1311,1415,1901,1906,1911,503,1312,1413,515)
outlierTrue <- sort(unique(outlierTrue))
label.true <- rep(0, wafer.amount)
label.true[wafer.id%in%outlierTrue] <- 1
#### Use die mean for the analysis
ibm.data <- as.matrix(cbind(id, ibm[,-c(1:2)]))
die.avg <- rowMeans(ibm.data[,4:ncol(ibm.data)])
die.sd <- apply(ibm.data[,4:ncol(ibm.data)], 1, sd)
plot(die.avg, die.sd, pch=19)
wafer.avg <- ave(die.avg, as.factor(ibm.data[,1]))
w.d.res <- matrix(die.avg-wafer.avg, wafer.amount, byrow=T, dimnames=list(wafer.id, die.loc))
# Each element in w.d.res represents the mean of 14 measures on a given die.
# Each row and column represents a wafer and a certain position.
die.pca <- prcomp(ibm.data[,-(1:3)], T, T, F) # PCA on die level
cumsum(die.pca$sdev^2)/sum(die.pca$sdev^2) # First component takes over 85% of the variance
die.pca$rot[,1] # Coefficients are similar for each measurement, so using mean is valid
#### Divide data into training and test sets
train.size <- round(wafer.amount*0.8)
test.size <- wafer.amount-train.size
train.id <- sample(wafer.id, train.size)
train.id <- train.id[order(train.id)]
test.id <- wafer.id[wafer.id %in% train.id == F]
ibm.train <- ibm.data[id %in% train.id,]
ibm.test <- ibm.data[id %in% test.id,]
#### PCA on the die residuals
pca.r <- prcomp(w.d.res, T, T ,F)
cumsum(pca.r$sdev^2)/sum(pca.r$sdev^2) # First 3 components contain 62% of the total variance
label.pca.r <- rep(F, wafer.amount)
label.pca.r[pca.r$x[,1]>3 | pca.r$x[,1]<(-3)] <- T
label.pca.r[pca.r$x[,2]>2.5 | pca.r$x[,2]<(-1.5)] <- T
label.pca.r[pca.r$x[,3]>2 | pca.r$x[,3]<(-2)] <- T
png("Figure1.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(pca.r$x[,1], pca.r$x[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component", col=label.pca.l+1)
abline(h=c(-1.5, 2.5), v=c(-3, 3), lwd=2, col=c(4,4,2,2))
plot(pca.r$x[,2], pca.r$x[,3], pch=19, xlab="Second Principal Component", ylab="Third Principal Component", col=label.pca.l+1)
abline(h=c(-2, 2), v=c(-1.5, 2.5), lwd=2, col=c(6,6,4,4))
dev.off()
get.sum <- function(x){
  tot <- sum(x) # Number of fault wafers identified
  FR <- sum((wafer.id[x==T] %in% outlierTrue) == F) # False rate
  mis <- sum((outlierTrue %in% wafer.id[x==T]) == F) # Misclassification rate
  return(list(c("tot"=tot, "FR"=FR, "mis"=mis), wafer.id[x==T]))
}
get.sum(label.pca.r)
#### Emiprical variograms on complete dataset
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
#### PCA on variograms
pca.v <- prcomp(geo.v, T, T, F)
cumsum(pca.v$sdev^2)/sum(pca.v$sdev^2) # First 3 components contain 98% of the total variance
label.pca.v <- rep(F, wafer.amount)
label.pca.v[pca.v$x[,1]>1] <- T
label.pca.v[pca.v$x[,2]>0.2 | pca.v$x[,2]<(-0.4)] <- T
label.pca.v[pca.v$x[,3]>0.15 | pca.v$x[,3]<(-0.2)] <- T
png("Figure2.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(pca.v$x[,1], pca.v$x[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component", col=label.pca.v+1)
abline(h=c(-0.4, 0.2), v=1, lwd=2, col=c(4,4,2))
plot(pca.v$x[,2], pca.v$x[,3], pch=19, xlab="Second Principal Component", ylab="Third Principal Component", col=label.pca.v+1)
abline(h=c(-0.2, 0.15), v=c(-0.4, 0.2), lwd=2, col=c(6,6,4,4))
dev.off()
get.sum(label.pca.v)
#### Die residuals, normal kernel PCA
kpca.r <- kpca(w.d.res, kpar=list(sigma=0.5))
rot.k.r <- rotated(kpca.r)
rot.k.r2 <- rotated(kpca(w.d.res, kpar=list(sigma=1)))
label.kpca.r <- label.kpca.r2 <- rep(F, wafer.amount)
label.kpca.r[rot.k.r[,1]>1.5 & rot.k.r[,2]<1.5] <- T
label.kpca.r2[rot.k.r2[,1]<(-0.5) & rot.k.r2[,2]>(-0.5)] <- T
png("Figure3.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(rot.k.r[,1], rot.k.r[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component", col=label.kpca.r+1)
plot(rot.k.r2[,1], rot.k.r2[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component", col=label.kpca.r2+1)
dev.off()
get.sum(label.kpca.r); get.sum(label.kpca.r2)
#### Variogram, normal kernel PCA
kpca.v <- kpca(geo.v, kpar=list(sigma=1))
rot.k.v <- rotated(kpca.v)
label.kpca.v <- rep(F, wafer.amount)
label.kpca.v[rot.k.v[,1]>10] <- T
label.kpca.v[rot.k.v[,2]>3] <- T
png("Figure4.png", 1600, 800)
par(mfrow=c(1,2), mar=c(5,5,2,3), cex=1.8, cex.lab=1.5, cex.axis=1.3)
plot(rot.k.v[,1], rot.k.v[,2], pch=19, xlab="First Principal Component", ylab="Second Principal Component",  col=label.kpca.v+1)
plot(rot.k.v[,2], rot.k.v[,3], pch=19, xlab="Second Principal Component", ylab="Third Principal Component",  col=label.kpca.v+1)
dev.off()
get.sum(label.kpca.v)
