rm(list=ls())
require(geoR)
load('ibm.Rdata')

wid <- unique(id)
dieTT <- table(id)
notNA <- as.numeric(names(dieTT)[which(dieTT == max(dieTT))])
nna <- lw.id(notNA[2])
nna.x <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipX']
nna.y <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipY']
i=1;j=1;

meas <- grep("psro", names(ibm))
sill <- matrix(0, nrow=length(unique(id)), ncol=length(meas))
#for (i in 1:length(wid)){
#  for (j in 1:length(meas)){
    gid <- lw.id(wid[i])
    gR <- ibm[ibm$lot == gid$l & ibm$wafer == gid$w,  ]
    gR <- as.geodata(gR, coords.col=c('chipX','chipY'), data.col=meas[j])
    vario <- variog(gR, estimate='modulus')
    var.mod <- variofit(vario, c(max(vario$v),max(vario$u)), cov.model='cubic')
    sill[i, j] <- var.mod$cov.pars[1]
#  }
#}
load('sill.Rdata')
cid <- which(sill == 0, arr=T)
table(cid[  ,1])


sill2 <- rep(0, length(wid))
#for (i in 1:length(wid)){
  gid <- lw.id(wid[i])
  gR <- ibm[ibm$lot == gid$l & ibm$wafer == gid$w,  ]
  dieMeans <- rowMeans(gR[  ,meas])
  gR <- cbind(gR,dieMeans)
  gR <- as.geodata(gR, coords.col=c('chipX','chipY'), data.col='dieMeans')
  vario <- variog(gR, estimate='modulus')
  var.mod <- variofit(vario, c(max(vario$v),max(vario$u)), cov.model='cubic')
  sill2[i] <- var.mod$cov.pars[1]

# Grphical exam
#  fn <- sprintf('Lot%dWafer%d.jpeg',gid$l,gid$w)
#  jpeg(file=fn)
#  plot(vario,main='Variogram');lines(var.mod)
#  dev.off()
#}
load('sill2.Rdata')
table(cid[cid[ ,1] %in% which(sill2 == 0), 1])
which(sill2 == 0)


sill.est <- sill
for (k in 1:dim(cid)[1]){
  if (sill2[cid[k,1]] != 0)
    sill.est[cid[k,1],cid[k,2]] <- sill2[cid[k,1]]
  else{
    ss <- sill[cid[k,1],  ]
    ss <- ss[ss != 0]
    sill.est[cid[k,1],cid[k,2]] <- median(ss)
  }
}
which(sill.est == 0)






wfgr <- data.frame(cbind(nna.x, nna.y))
names(wfgr) <- c('chipX','chipY')
rng <- min(max(nna.x),max(nna.y))
nug <- 0

ndie <- length(nna.x)
IBMC <- data.frame(matrix(0, nrow=length(wid)*ndie,
                          ncol=ncol(ibm)))
names(IBMC) <- names(ibm)

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


## wafer meas heatmap
plotWmap <- function(dm, id){
  title <- sprintf("Lot%dWafer%d", id$l, id$w)
  fig <- ggplot(dm, aes(x = chipX, y = chipY)) + scale_x_continuous(breaks=1:13) + scale_y_continuous(breaks=1:13)
  fig <- fig + geom_point(aes(colour = intpol),shape = 15,size = 10) + ggtitle(title)
  fig <- fig + guides(colour = guide_colorbar(barwidth = 0.5, barheight = 20))
  fig
}
dm <- data.frame(cbind(wfgr,intpol))
plotWmap(dm, gid)

