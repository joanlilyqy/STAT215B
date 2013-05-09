#### STAT215B Project
#### Ying Qiao
require(ggplot2)
require(geoR)


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

