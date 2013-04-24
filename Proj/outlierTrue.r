rm(list=ls())
load('ibm.Rdata')

wid <- unique(id)
dieTT <- table(id)
notNA <- as.numeric(names(dieTT)[which(dieTT == max(dieTT))])
nna <- lw.id(notNA[2])
nna.x <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipX']
nna.y <- ibm[ibm$lot == nna$l & ibm$wafer == nna$w, 'chipY']
i=1;j=1;k=1


## wafer meas heatmap
plotWmap <- function(dm, id){
  title <- sprintf("Lot%dWafer%d", id$l, id$w)
  fn <- sprintf("%s.jpeg",title)
  fig <- ggplot(dm, aes(x = chipX, y = chipY)) + scale_x_continuous(breaks=1:13) + scale_y_continuous(breaks=1:13)
  fig <- fig + geom_point(aes(colour = dieMean),shape = 15,size = 10) + ggtitle(title)
  fig <- fig + guides(colour = guide_colorbar(barwidth = 0.5, barheight = 20))
  fig
#  ggsave(filename=fn, plot=fig)
}

meas <- grep("psro", names(ibm))
wfgr <- data.frame(cbind(nna.x, nna.y))
names(wfgr) <- c('chipX','chipY')
medianww <- data.frame(matrix(0, nrow=length(nna.x), ncol=length(meas)))
names(medianww) <- names(ibm)[meas]
for (k in 1:length(nna.x)){
  xx <- wfgr[k, 'chipX']
  yy <- wfgr[k, 'chipY']
  wwData <- ibm[ibm$chipX == xx & ibm$chipY == yy,  ]
  medianww[k,  ] <- apply(wwData[  , meas], 2, median)
}
goldwafer <- data.frame(cbind(wfgr, medianww))

dm <- goldwafer[ ,1:3]
fn <- "gold_wafer.jpg"
fig <- ggplot(dm, aes(x = chipX, y = chipY)) + scale_x_continuous(breaks=1:13) + scale_y_continuous(breaks=1:13)
fig <- fig + geom_point(aes(colour = dm[ ,3]),shape = 15,size = 10) + ggtitle("Gold Median Wafer")
fig <- fig + guides(colour = guide_colorbar(barwidth = 0.5, barheight = 20))
#fig
#ggsave(filename=fn, plot=fig)

    
gR <- as.geodata(dm, coords.col=c('chipX','chipY'), data.col=3)
vario <- variog(gR, estimate='modulus')
var.mod <- variofit(vario, c(0.1, 6), cov.model='cubic')
#jpeg(file='vario_fit.jpg')
plot(vario,main='Empirical Variogram and Estimated Model');lines(var.mod)
#dev.off()





#for (i in wid){
  pid <- lw.id(i)
  wwData <- ibm[ibm$lot == pid$l & ibm$wafer == pid$w,  ]
  dieMean<-rowMeans(wwData[  ,grep("psro", names(ibm))])
  dm <- cbind(wwData[  ,c('chipX','chipY')], dieMean)
#  plotWmap(dm, pid)
#}


outlierTrue <- c(421,507,523,707,1605,1607,2201,2205,508,703,421,
                 1305,1412,1507,1604,1606,1207,1201,607,514,511,509,505,1909,1915,
                 1705,517,520,521,1217,1311,1415,1901,1906,1911,503,1312,1413,515)
outlierTrue <- sort(unique(outlierTrue))
length(outlierTrue)

