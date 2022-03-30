
library(dtwclust)
library(dplR)
library(TSdist)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import combined fluorescent data 
Data = read.csv('AllDeidentified081921.csv')  
Data = Data[1:4799,]
Data = Data[-c(1)]
DataT = t(Data)


#####################
# Data preprocess
#####################

# ----------------------------
# Remove high frequency noise
# ----------------------------
Data.LowP = apply(DataT, 1, 
                  function(x)pass.filt(x, W=0.05, type="low", method="Butterworth"))
plot(Data$X0682.0373[1:3600], type="l", col="grey10",ylim=c(-5,45))#,axes=FALSE) 
abline(h=0, col="grey70",lty="dashed",lwd=2.0)
lines(Data.LowP[1:3600,"X0682.0373"], type="l", col = "cornflowerblue", lwd=3.0)

# ----------------------------
# Remove baseline drift 
# ----------------------------
Data.NoInd <- Data.LowP
Data.trend <- Data.LowP
Data.NoInd[c(500:750,1200:1350,1800:1950,2400:2550,3000:3150,3600:3800,4200:4400),] = NaN 
for (i in colnames(Data.LowP)){
  for (j in 1:8){ 
    if (j == 1){
      Data.NoInd[1:(j*600),i] = pass.filt(t(Data.LowP[1:(j*600),i]), W=0.005, type="low", method="Butterworth")
    }
    if (j>1 & j<7){
      Data.NoInd[((j-1)*600+150):(j*600),i] = pass.filt(t(Data.LowP[((j-1)*600+150):(j*600),i]), W=0.005, type="low", method="Butterworth")
    }
    if (j==7){
      Data.NoInd[((j-1)*600+200):(j*600),i] = pass.filt(t(Data.LowP[((j-1)*600+200):(j*600),i]), W=0.005, type="low", method="Butterworth")
    }
    if (j==8){
      Data.NoInd[((j-1)*600+200):(nrow(Data)),i] = pass.filt(t(Data.LowP[((j-1)*600+200):(nrow(Data)),i]), W=0.005, type="low", method="Butterworth")
    }
  }
  Data.trend[,i] = spline(Data.NoInd[,i], method = "natural", n=nrow(Data))$y   
}
plot(Data.LowP[1:3600,"X0682.0373"], type="l", col="grey10", ylim=c(-5,45))#,axes=FALSE)
abline(h=0, col="grey70", lty="dashed", lwd=2.0)
lines(Data.NoInd[1:3600,"X0682.0373"], type="l", col="indianred2", lwd=5.0) 
lines(Data.trend[1:3600,"X0682.0373"], col="cornflowerblue", lwd=3.0)  

Data.0Mean <- Data.LowP-Data.trend
DataT.0Mean = t(Data.0Mean)
plot(Data.0Mean[1:3600,"X0682.0373"],type="l", col="grey10",ylim=c(-5,45))#,axes=FALSE)
abline(h=0, col="grey70",lty="dashed", lwd=2.0)


# ----------------------------
# Response detection
# ----------------------------
Data.Rsp <- Data.0Mean
basetime <- 580
baseparts <- 1
baseinter <- basetime/baseparts
for (i in 1:ncol(Data.0Mean)){
  FreFlted <- rep(0, nrow(Data.0Mean))
  maxsum <- 0.0
  for (ibase in 1:baseparts){
    tempmax <- max(Data.0Mean[(ibase*baseinter-baseinter+1):(ibase*baseinter),i])
    maxsum <- maxsum + tempmax
  }
  QHigh <- maxsum*2/baseparts 
  AmpFlted = FreFlted
  for (j in 1:nrow(Data.0Mean)){
    if (Data.0Mean[j,i] > QHigh){
      AmpFlted[j] = Data.0Mean[j,i] - QHigh 
    }
  }
  Data.Rsp[,i] <- AmpFlted
}
plot(Data.0Mean[1:3600,"X0682.0373"], type="l", col="grey10",ylim=c(-5,45))#,axes=FALSE)  
abline(h=0, col="grey70", lty="dashed", lwd=2.0) 
lines(Data.Rsp[1:3600,"X0682.0373"], col="cornflowerblue", lwd=2.0)
plot(Data.Rsp[1:3600,"X0682.0373"], type="l", col="grey10",ylim=c(-5,45))#,axes=FALSE)  
abline(h=0, col="grey70", lty="dashed", lwd=2.0) 

# ----------------------------
# Normalization
# ----------------------------
DataN <- Data.Rsp[1:3600,]  # only responses to indentation stimuli (brushing excluded)
MaxForAll <- max(DataN)
for (i in colnames(DataN)){
  DataN[,i] = Data.Rsp[1:3600,i]*max(Data.Rsp[1:3600,])/max(Data.Rsp[1:3600,i])
}
DataN[is.na(DataN)] <- 0
DataNT = t(DataN)
plot(Data.Rsp[1:3600,"X0682.0373"], type="l", col="firebrick", ylim=range(-5:45), lwd=2)
lines(DataN[,"X0682.0373"], type="l", col="dodgerblue3")
plot(DataN[,"X0682.0373"], type="l", col="grey10",ylim=range(-5:45))
abline(h=0, col="grey70", lty="dashed", lwd=2) 

# ----------------------------
# Remove long intervals
# ----------------------------
Inden = DataN[1:3600,]
IndenT <- t(Inden)
Inden.NoIntv <- DataN[c(500:800,1200:1400,1800:2000,2400:2600,3000:3200),]
IndenT.NoIntv = t(Inden.NoIntv)
plot(Inden.NoIntv[,"X0682.0373"], type="l", col="grey10",ylim=range(-5:45))  
abline(h=0, col="grey70", lty="dashed", lwd=2) 



########################
# Hierarchical clustering
########################

# register distance metric: FourierDistance
proxy::pr_DB$set_entry(FUN = FourierDistance, names = c("Fourier"),
                       loop = TRUE, type = "metric", distance = TRUE,
                       description = "Customized distance")

# ----------------------------
# 1st layer
# ----------------------------
hc <- tsclust(IndenT.NoIntv, type = "hierarchical",
              k =4, trace = TRUE, distance = "Fourier", 
              control = hierarchical_control(method = "ward.D"))
plot(hc, type = "series")
hc@clusinfo
#hc_idx <- hc@cluster; names(hc_idx) <- NULL; hc_idx

# ----------------------------
# 2nd layer
# ----------------------------
# data before response detection
Inden = Data.0Mean[1:3600,]
IndenT <- t(Inden)
Inden.NoIntv <- Data.0Mean[c(500:800,1200:1400,1800:2000,2400:2600,3000:3200),]  
IndenT.NoIntv = t(Inden.NoIntv)

group1_idx <- c(which(hc@cluster==1))
IndenT.NoIntv.group1 <- IndenT.NoIntv[group1_idx,]
hc.group1 <- tsclust(IndenT.NoIntv.group1, type = "hierarchical",
                    k = 4, trace = TRUE, distance = "Fourier",
                    control = hierarchical_control(method = "ward.D"))
plot(hc.group1, type = "sc")
hc.group1@clusinfo
#hc.group1_idx <- hc.group1@cluster; names(hc.group1_idx) <- NULL; hc.group1_idx

group2_idx <- c(which(hc@cluster==2))
IndenT.NoIntv.group2 <- IndenT.NoIntv[group2_idx,]
hc.group2 <- tsclust(IndenT.NoIntv.group2, type = "hierarchical",
                    k = 4, trace = TRUE, distance = "Fourier",
                    control = hierarchical_control(method = "ward.D"))
plot(hc.group2, type = "sc")
hc.group2@clusinfo
#hc.group2_idx <- hc.group2@cluster; names(hc.group2_idx) <- NULL; hc.group2_idx

group3_idx <- c(which(hc@cluster==3))
IndenT.NoIntv.group3 <- IndenT.NoIntv[group3_idx,]
hc.group3 <- tsclust(IndenT.NoIntv.group3, type = "hierarchical",
                     k = 4, trace = TRUE, distance = "Fourier",
                     control = hierarchical_control(method = "ward.D"))
plot(hc.group3, type = "sc")
hc.group3@clusinfo
#hc.group3_idx <- hc.group3@cluster; names(hc.group3_idx) <- NULL; hc.group3_idx

group4_idx <- c(which(hc@cluster==4))
IndenT.NoIntv.group4 <- IndenT.NoIntv[group4_idx,]
hc.group4 <- tsclust(IndenT.NoIntv.group4, type = "hierarchical",
                     k = 4, trace = TRUE, distance = "Fourier",
                     control = hierarchical_control(method = "ward.D"))
plot(hc.group4, type = "sc")
hc.group4@clusinfo
#hc.group4_idx <- hc.group4@cluster; names(hc.group4_idx) <- NULL; hc.group4_idx

