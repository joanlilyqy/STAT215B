rm(list=ls())
load('ibm.Rdata')


## wafer meas heatmap
plotWmap <- function(dm, id){
  title <- sprintf("Lot%dWafer%d", id$l, id$w)
  fn <- sprintf("%s.jpeg",title)
  fig <- ggplot(dm, aes(x = chipX, y = chipY)) + scale_x_continuous(breaks=1:13) + scale_y_continuous(breaks=1:13)
  fig <- fig + geom_point(aes(colour = dieSD),shape = 15,size = 10) + ggtitle(title)
  fig <- fig + guides(colour = guide_colorbar(barwidth = 0.5, barheight = 20))
  fig
#  ggsave(filename=fn, plot=fig)
}


ibm$idno <- ibm$lot_ID*100 + ibm$wafer_ID
for (mi in unique(ibm$idno)){
mi=707
nna <- lw.id(mi)
plotid <- nna
wwData <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w,  ]
dieMean<-rowMeans(wwData[  ,grep("psro", names(ibm))])
dieSD <-rowSums((wwData[  ,grep("psro", names(ibm))] - dieMean)^2)/14
dm <- cbind(wwData[  ,c('chipX','chipY')], dieSD) #dieMean, dieSD
plotWmap(dm, nna)
#}






ll <- c(707,523,2205,301,1008,801,1005,2201,1607,914,
        1619,1620,2316,807,808,911,1512,1513,1605,1617,1618,2306,2307,2311,2314,2315)
