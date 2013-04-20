#### STAT215B Project ####

rm(list=ls())
ibm <- read.table("IBMBigData.csv", sep=",", header=T)
dim(ibm)
n <- dim(ibm)[2] - 4
head(ibm)

require(ggplot2)
## Fig.1
table(ibm$lot_ID)
id <- ibm$lot_ID*100 + ibm$wafer_ID
length(unique(id))
sum(table(floor(unique(id)/100)))

plot(table(floor(unique(id)/100)),col='blue',lwd=2,typ='h',xlab='Lot ID',ylab='Num. of wafers within lot')

## Fig.2
l23w16 <- ibm[ibm$lot_ID == 23 & ibm$wafer_ID == 16,]
dim(l23w16)
dieMean<-rowMeans(l23w16[,5:18])
dm <- cbind(l23w16[,3:4],dieMean)

fig <- ggplot(dm, aes(x = chipX, y = chipY)) + scale_x_continuous(breaks=1:13) + scale_y_continuous(breaks=1:13)
fig <- fig + geom_point(aes(colour = dieMean),shape = 15,size = 10) + ggtitle("Lot 23 | Wafer 16")
fig + guides(colour = guide_colorbar(barwidth = 0.5, barheight = 20))

