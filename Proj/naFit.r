rm(list=ls())
require(geoR)
load('ibm.Rdata')

wid <- unique(id)
dieTT <- table(id)
notNA <- as.numeric(names(dieTT)[which(dieTT == max(dieTT))])
nna <- lw.id(notNA[2])
nna.x <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipX']
nna.y <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipY']

meas <- grep("psro", names(ibm))
sill <- matrix(0, nrow=length(unique(id)), ncol=length(meas))
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
cid <- which(sill == 0, arr=T)
table(cid[  ,1])

sill2 <- rep(0, length(wid))
for (i in 1:length(wid)){
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
}
wid[which(sill2 == 0)]

sort(table(cid[cid[ ,1] %in% which(sill2 == 0), 1]))
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
sill2.est <- sill2
sill2.est[sill2.est == 0] <- apply(sill.est[which(sill2 == 0), ], 1, median)




