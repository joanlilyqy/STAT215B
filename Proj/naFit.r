#rm(list=ls())
#load('ibm.Rdata')

dieTT <- table(id)
notNA <- as.numeric(names(dieTT)[which(dieTT == max(dieTT))])
nna <- lw.id(notNA[1])
x <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipX']
y <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipY']

dd <- nna # ?id
ww <- ibm[ibm$lot == dd$l & ibm$wafer == dd$w, 'psro_00']
cx <- 6#x[which(ww == max(ww))]
cy <- 6#y[which(ww == max(ww))]
x1 <- x - cx
y1 <- y - cy
XY <- cbind(poly(x1, 6, raw=TRUE), poly(y1, 6, raw=TRUE), poly(x1*y1, 3, raw=TRUE))
fit <- lm(ww ~ XY)
summary(fit)
summary(fit)$sigma

## wafer meas heatmap
plotWmap <- function(dm, id){
  title <- sprintf("Lot %d | Wafer %d", id$l, id$w)
  fig <- ggplot(dm, aes(x = chipX, y = chipY)) + scale_x_continuous(breaks=1:13) + scale_y_continuous(breaks=1:13)
  fig <- fig + geom_point(aes(colour = dm[  ,3]),shape = 15,size = 10) + ggtitle(title)
  fig + guides(colour = guide_colorbar(barwidth = 0.5, barheight = 20))
}

plotid <- nna
wwData <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w,  ]
dieMean<-rowMeans(wwData[  ,grep("psro", names(ibm))])
dm <- cbind(wwData[  ,c('chipX','chipY')], fit$res)
plotWmap(dm, nna)

