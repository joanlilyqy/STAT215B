rm(list=ls())
load('ibm.Rdata')


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

wid <- unique(id)
for (mi in wid){
  pid <- lw.id(mi)
  wwData <- ibm[ibm$lot == pid$l & ibm$wafer == pid$w,  ]
  dieMean<-rowMeans(wwData[  ,grep("psro", names(ibm))])
  dm <- cbind(wwData[  ,c('chipX','chipY')], dieMean)
  plotWmap(dm, pid)
}


outlierTrue <- c(421,507,523,707,1605,1607,2201,2205,508,703,421,
                 1305,1412,1507,1604,1606,1207,1201,607,514,511,509,505,1909,1915,
                 1705,517,520,521,1217,1311,1415,1901,1906,1911,503,1312,1413,515)
outlierTrue <- sort(unique(outlierTrue))
length(outlierTrue)

