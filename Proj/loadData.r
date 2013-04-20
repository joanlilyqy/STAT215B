#### STAT215B Project ####

rm(list=ls())
require(ggplot2)


ibm <- read.table("IBMBigData.csv", sep=",", header=T)
dim(ibm)
head(ibm)
table(ibm$lot_ID)


id <- ibm$lot_ID*100 + ibm$wafer_ID
M <- length(unique(id))
sum(table(floor(unique(id)/100)))


lw.id <- function(ii){
  lid <- floor(ii/100)
  wid <- ii - lid*100
  return (list(lid=lid,wid=wid))
}

save(list=ls(), file='ibm.Rdata')

