# Haotian's Code
rm(list=ls())
setwd("~/Stat215B/Final Project")
library(geoR)
ibm <- read.table("IBMC.csv", sep=",", header=T)
library(e1071)

ibm <- ibm[,-1]
id <- ibm$lot_ID*100 + ibm$wafer_ID # Define a new id for each wafer
loc <- ibm$chipX*100 + ibm$chipY # Define new positions for each die
wafer.amount <- length(unique(id)) # Total amount of wafers
wafers <- unique(id)

# average 14 observations
ibm.MeanChips<- ibm[,1:4]
ibm.MeanChips$loc <- ibm.MeanChips$chipX*100 + ibm.MeanChips$chipY
ibm.MeanChips$psromean <- rowMeans(ibm[,5:18])
ibm.NoSpatial <- matrix(0,wafer.amount,117+2)
for (iWafer in 1:wafer.amount)
{
  data <- ibm.MeanChips[id==wafers[iWafer],]
  Order <- order(data$loc)
  ibm.NoSpatial[iWafer,]<-c(floor(wafers[iWafer]/100),
    wafers[iWafer]-floor(wafers[iWafer]/100)*100,data[Order,6])
}
rm(data)
rm(Order)
colnames(ibm.NoSpatial)<-c("lot_ID","wafer_ID",as.character(sort(unique(loc))))
ibm.NoSpatial = as.data.frame(ibm.NoSpatial)

# outliers
outlierTrue <- c(421,507,523,707,1605,1607,2201,2205,508,703,421,
                 1305,1412,1507,1604,1606,1207,1201,607,514,511,
                 509,505,1909,1915,1705,517,520,521,1217,1311,
                 1415,1901,1906,1911,503,1312,1413,515)
outlierTrue <- sort(unique(outlierTrue))
outlierTrue.Matrix = cbind(floor(outlierTrue/100),
                           outlierTrue-floor(outlierTrue/100)*100)



# treat response as a binary factor, with 1 equal to yes outlier and 0 no.
ibm.outlier <- apply(ibm.NoSpatial[,1:2],1,function(x){sum(colSums(x==t(outlierTrue.Matrix))==2)})
Fault <- c("Normal","Fault")
ibm.NoSpatial$Fault <- as.factor(Fault[ibm.outlier+1])
head(ibm.NoSpatial)

# create the one with spatial correlation
ibm.Spatial <- matrix(0,nrow=wafer.amount,ncol=2+12)
for (iWafer in 1:wafer.amount)
{
  data <- ibm.MeanChips[id==wafers[iWafer],]
  data <- as.geodata(data, coords.col=c('chipX','chipY'),data.col=6)
  ibm.Spatial[iWafer,]=c(floor(wafers[iWafer]/100),
    wafers[iWafer]-floor(wafers[iWafer]/100)*100,
    variog(data, estimate='modulus',uvec=seq(0,12,by=1))$v)
}
rm(data)
rm(iWafer)

ibm.Spatial <- as.data.frame(ibm.Spatial)
ibm.outlier <- apply(ibm.Spatial[,1:2],1,function(x){sum(colSums(x==t(outlierTrue.Matrix))==2)})
ibm.Spatial$Fault <- as.factor(Fault[ibm.outlier+1])
rm(Fault)

# train and test
# a function that calculates the predictions by taking in data, a vector that
# specifies which data points are for training, the parameters and the method
Pars.SVM <- function(data,train,pars,method)
{
  # cost and gamma are two vectors
  results <- apply(pars,1,function(x){Model<-svm(Fault~.,data=data[train,-(1:2)],
                       type=method,cost=x[1],gamma=x[2]);
  as.numeric(table(predict(Model,data[-train,-c(1,2,dim(data)[2])]),
        data[-train,dim(data)[2]]));})
  # colnames(results)<-1:dim(pars)[1]
  # rownames(results)<-c("TP","FN","FP","TN")
  results = t(results)
  # TP <- which(results[,1]==max(results[,1]))
  # TPFP <- which(results[TP,3]==min(results[TP,3]))
  # if (length(TPFP)==1)
  # {results = c(results[TP[TPFP],],TP[TPFP])}
  # else
  # {results = cbind(results[TP[TPFP],],TP[TPFP])}
  return(results)
}

cost = 10^(seq(1,4,length.out=30))
gamma = seq(0.001,0.030,by=0.001)
pars = matrix(matrix(1,nrow=length(gamma),ncol=length(cost))%*%diag(cost),ncol=1)
pars = cbind(pars, matrix(diag(gamma)%*%matrix(1,nrow=length(gamma),ncol=length(cost)),ncol=1))

# single dies and c
# C.Class = matrix(0,nrow=100,ncol=4)
# for (momo in 100:100)
# {
#   train <- sample(1:wafer.amount,floor(wafer.amount*0.5))
#   C.Class[momo,] = Pars.SVM(ibm.NoSpatial,train,pars,"C-classification")[1,1:4]
# }
train <- sample(1:wafer.amount,floor(wafer.amount*0.5))
Die.C = Pars.SVM(ibm.NoSpatial,train,pars,"C-classification")
Die.C.Tradeoff = cbind(Die.C[,1]/(Die.C[,1]+Die.C[,2]),Die.C[,1]/(Die.C[,1]+Die.C[,3]))
Die.C.Tradeoff = Die.C.Tradeoff[order(Die.C.Tradeoff[,1]),]
# write.table(C.Class,"C-classification.csv",sep=",",row.names=F,col.names=F)
# par(mfrow=c(2,1))
# hist(C.Class[,1]/rowSums(C.Class[,1:2]),main="Rate of Correct Prediction among All Faults - Single Dies",
#      xlab="Correct Prediction Rate",breaks=seq(0,1,by=0.05))
# hist(C.Class[,3]/rowSums(C.Class[,3:4]),main="Rate of False Prediction among All Normals - Single Dies",
#      xlab="False Prediction Rate",breaks=seq(0,1,by=0.05))

# vario and c
# C.Class.Spatial = matrix(0,nrow=100,ncol=4)
# for (momo in 1:100)
# {
#   train <- sample(1:wafer.amount,floor(wafer.amount*0.5))
#   C.Class.Spatial[momo,] = Pars.SVM(ibm.Spatial,train,pars,"C-classification")[1,1:4]
# }
Vario.C = Pars.SVM(ibm.Spatial,train,pars,"C-classification")
Vario.C.Tradeoff = cbind(Vario.C[,1]/(Vario.C[,1]+Vario.C[,2]),Vario.C[,1]/(Vario.C[,1]+Vario.C[,3]))
Vario.C.Tradeoff = Vario.C.Tradeoff[order(Vario.C.Tradeoff[,1]),]
# write.table(C.Class.Spatial,"C-classification-Spatial.csv",sep=",",row.names=F,col.names=F)
# par(mfrow=c(2,1))
# hist(C.Class.Spatial[,1]/rowSums(C.Class.Spatial[,1:2]),main="Rate of Correct Prediction among All Faults - Variogram",
#      xlab="Correct Prediction Rate",breaks=seq(0,1,by=0.05))
# hist(C.Class.Spatial[,3]/rowSums(C.Class.Spatial[,3:4]),main="Rate of False Prediction among All Normals - Variogram",
#      xlab="False Prediction Rate",breaks=seq(0,1,by=0.05))
# 

Novelty.SVM <- function(data,train,pars,Method)
{
  results <- apply(pars,1,function(x){Model<-svm(Fault~.,data=data[train,-(1:2)],
                                                 type=Method,gamma=x[2],cost=x[1]);
                                      as.numeric(table(predict(Model,data[-train,-c(1,2,dim(data)[2])]),
                                                       data[-train,dim(data)[2]]));})
  # colnames(results)<-1:dim(pars)[1]
  # rownames(results)<-c("TP","FN","FP","TN")
  results = t(results)
#   TP <- which(results[,1]==max(results[,1]))
#   TPFP <- which(results[TP,3]==min(results[TP,3]))
#   if (length(TPFP)==1)
#     {results = c(results[TP[TPFP],],TP[TPFP])}
#   else
#     {results = cbind(results[TP[TPFP],],TP[TPFP])}
  return(results)
}


Die.Novelty <- Novelty.SVM(ibm.NoSpatial,train,pars,"one-classification")
Die.Novelty.TO <- cbind(Die.Novelty[,1]/(Die.Novelty[,1]+Die.Novelty[,2]),Die.Novelty[,1]/(Die.Novelty[,1]+Die.Novelty[,3]))
Die.Novelty.TO <- Die.Novelty.TO[order(Die.Novelty.TO[,1]),]

Vario.Novelty <- Novelty.SVM(ibm.Spatial,train,pars,"one-classification")
Vario.Novelty.TO <- cbind(Vario.Novelty[,1]/(Vario.Novelty[,1]+Vario.Novelty[,2]),Vario.Novelty[,1]/(Vario.Novelty[,1]+Vario.Novelty[,3]))
Vario.Novelty.TO <- Vario.Novelty.TO[order(Vario.Novelty.TO[,1]),]

plot(Die.C.Tradeoff[,2],Die.C.Tradeoff[,1],xlab="Precision",ylab='Recall',
     main="Comparison of Trade-off between Precision and Recall",
     type="l",xlim=c(0.15,1),ylim=c(0,1),cex.lab=1.5,cex.main=1.5)

lines(Vario.C.Tradeoff[,2],Vario.C.Tradeoff[,1],lty=2,lwd=2.5)
lines(Vario.Novelty.TO[,2],Vario.Novelty.TO[,1],lty=2,lwd=2.5)
lines(Die.Novelty.TO[,2],Die.Novelty.TO[,1],lty=2)
# library(plotrix)
draw.circle(mean(range(Vario.Novelty.TO[,2])),
            mean(Vario.Novelty.TO[,1]),
            diff(range(Vario.Novelty.TO[,2]))/2)
legend('topright',c('Die+C','Vario+C','Die+1','Vario+1'),
       lty=c(1,2,2,2),lwd=c(1,2.5,1,2.5))
# # variogram and one-class
# OneClass.Spatial <- matrix(0,nrow=100,ncol=4)
# for (momo in 1:100)
# {
#   train <- sample(1:wafer.amount,floor(wafer.amount*0.5))
#   temp<-Novelty.SVM(ibm.Spatial,train,gamma,"one-classification")
#   if (is.null(dim(temp)))
#   {OneClass.Spatial[momo,] = temp[1:4]}
#   else
#   {OneClass.Spatial[momo,] = temp[1,1:4]}
# }
# write.table(OneClass.Spatial,"1-classification-Spatial.csv",sep=",",row.names=F,col.names=F)
# par(mfrow=c(2,1))
# hist(OneClass.Spatial[,1]/rowSums(OneClass.Spatial[,1:2]),main="Rate of Correct Prediction among All Faults - Variogram",
#      xlab="Correct Prediction Rate",breaks=seq(0,1,by=0.05))
# hist(OneClass.Spatial[,3]/rowSums(OneClass.Spatial[,3:4]),main="Rate of False Prediction among All Normals - Variogram",
#      xlab="False Prediction Rate",breaks=seq(0,1,by=0.05))

# variogram and one-class
# OneClass <- matrix(0,nrow=100,ncol=4)
# for (momo in 1:100)
# {
#   train <- sample(1:wafer.amount,floor(wafer.amount*0.5))
#   temp<-Novelty.SVM(ibm.NoSpatial,train,gamma,"one-classification")
#   if (is.null(dim(temp)))
#   {OneClass[momo,] = temp[1:4]}
#   else
#   {OneClass[momo,] = temp[1,1:4]}
# }
# write.table(OneClass,"1-classification.csv",sep=",",row.names=F,col.names=F)
# par(mfrow=c(2,1))
# hist(OneClass[,1]/rowSums(OneClass[,1:2]),main="Rate of Correct Prediction among All Faults - Single Dies",
#      xlab="Correct Prediction Rate",breaks=seq(0,1,by=0.05))
# hist(OneClass[,3]/rowSums(OneClass[,3:4]),main="Rate of False Prediction among All Normals - Single Dies",
#      xlab="False Prediction Rate",breaks=seq(0,1,by=0.05))
# 
# boxplot(C.Class.Spatial[,1]/rowSums(C.Class.Spatial[,1:2]),C.Class[,1]/rowSums(C.Class[,1:2]),
#         OneClass.Spatial[,1]/rowSums(OneClass.Spatial[,1:2]),OneClass[,1]/rowSums(OneClass[,1:2]))
# boxplot(C.Class.Spatial[,3]/rowSums(C.Class.Spatial[,3:4]),C.Class[,3]/rowSums(C.Class[,3:4]),
#         OneClass.Spatial[,3]/rowSums(OneClass.Spatial[,3:4]),OneClass[,3]/rowSums(OneClass[,3:4]))
 
 
# C.Class.Both <- matrix(0,nrow=100,ncol=16)
# for (momo in 93:100)
# {
#   train <- sample(1:wafer.amount,floor(wafer.amount*0.5))
#   temp1<-Pars.SVM(ibm.NoSpatial,train,pars,"C-classification")
#   temp2<-Pars.SVM(ibm.Spatial,train,pars,"C-classification")
#   temp3<-Novelty.SVM(ibm.NoSpatial,train,gamma,"one-classification")
#   temp4<-Novelty.SVM(ibm.Spatial,train,gamma,"one-classification")
#   if (is.null(dim(temp1)))
#   {C.Class.Both[momo,1:4] = temp1[1:4]}
#   else
#   {C.Class.Both[momo,1:4] = temp1[1,1:4]}
#   if (is.null(dim(temp2)))
#   {C.Class.Both[momo,5:8] = temp2[1:4]}
#   else
#   {C.Class.Both[momo,5:8] = temp2[1,1:4]}
#   if (is.null(dim(temp3)))
#   {C.Class.Both[momo,9:12] = temp3[1:4]}
#   else
#   {C.Class.Both[momo,9:12] = temp3[1,1:4]}
#   if (is.null(dim(temp4)))
#   {C.Class.Both[momo,13:16] = temp4[1:4]}
#   else
#   {C.Class.Both[momo,13:16] = temp4[1,1:4]}
# }
# write.table(C.Class.Both,"all.csv",sep=",",row.names=F,col.names=F)
# par(mfrow=c(1,2))
# boxplot(C.Class.Both[,c(1,5,9,13)],ylab="Correction Predictions of Faults",
#         xlab="Method C=C-Classification, 1=novelty detection",
#         main="Comparison of Real Fault Detection",
#         names=c("Dies+C","Vario+C","Dies+1","Vario+1"))
# boxplot(C.Class.Both[,c(1,5,9,13)]/cbind(rowSums(C.Class.Both[,c(1:2)]),rowSums(C.Class.Both[,c(5:6)]),
#                                          rowSums(C.Class.Both[,c(9:10)]),rowSums(C.Class.Both[,c(13:14)])),
#         ylab="Rate Correction Predictions of Faults",
#         xlab="Method C=C-Classification, 1=novelty detection",
#         main="Comparison of Rates of Real Fault Detection",
#         names=c("Dies+C","Vario+C","Dies+1","Vario+1"))
# par(mfrow=c(1,2))
# boxplot(C.Class.Both[,c(3,7,11,15)],ylab="False Predictions of Normals",
#         xlab="Method C=C-Classification, 1=novelty detection",
#         main="Comparison of False Alarms",
#         names=c("Dies+C","Vario+C","Dies+1","Vario+1"))
# boxplot(C.Class.Both[,c(3,7,11,15)]/,/cbind(rowSums(C.Class.Both[,c(3:4)]),rowSums(C.Class.Both[,c(7:8)]),
#                                             rowSums(C.Class.Both[,c(11:12)]),rowSums(C.Class.Both[,c(15:16)])),
#         ylab="Rate of False Predictions of Normals",
#         xlab="Method C=C-Classification, 1=novelty detection",
#         main="Comparison of rates of False Alarms",
#         names=c("Dies+C","Vario+C","Dies+1","Vario+1"))
# t.test(C.Class.Both[,c(9,13)])
# t.test(C.Class.Both[,c(11,15)])