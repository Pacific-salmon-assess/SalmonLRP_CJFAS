

# library(tidyverse)
# library(ggplot2)
# library(viridis)
# library(corrplot)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
wcviCKDir<-paste(rootDir,"/WCVIChinookStudy",sep="")
setwd(wcviCKDir)

remove.EnhStocks <- TRUE
CoreInd <- FALSE#TRUE
AllExMH <- FALSE

# wcviCKSRDat <- read.csv(paste(wcviCKDir, "/DataIn/Inlet_Sum.csv", sep=""))
# wcviCKSRDat <- read.csv(paste(wcviCKDir, "/DataIn/Inlet_Sum_wEnh.csv", sep=""))

if(!CoreInd & !AllExMH){
  if(remove.EnhStocks) {wcviCKSRDat <- read.csv("DataIn/Inlet_Sum.csv")}
  if(!remove.EnhStocks) {wcviCKSRDat <- read.csv("DataIn/Inlet_Sum_wEnh.csv")}
}
if(CoreInd){wcviCKSRDat <- read.csv("DataIn/Inlet_Sum_CoreInd.csv")}
if(AllExMH){wcviCKSRDat <- read.csv("DataIn/Inlet_Sum_AllExMH.csv")}

wcviCKSRDat$yr_num <- group_by(wcviCKSRDat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
wcviCKSRDat$CU_ID <- group_by(wcviCKSRDat, Inlet_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
dum <- wcviCKSRDat %>% dplyr::select(CU_ID, BroodYear, Spawners)
# dum <- dum %>% tidyr::pivot_wider(id_cols=c(CU_ID, BroodYear), names_from=CU_ID,
#                            values_from=Spawners) %>% dplyr::select (!BroodYear)
dum <- dum %>% tidyr::pivot_wider(names_from=CU_ID,
                                  values_from=Spawners) %>% dplyr::select (!BroodYear)
dum <- dum %>% drop_na()
corMat <- cor(dum)

if(!CoreInd & !AllExMH){
  if(remove.EnhStocks) {
    rownames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars.csv",sep="/"))$stkName
    colnames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars.csv",sep="/"))$stkName
  }
  if(!remove.EnhStocks) {
    rownames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars_wEnh.csv",sep="/"))$stkName
    colnames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars_wEnh.csv",sep="/"))$stkName
  }
}
if(CoreInd){
  rownames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars_CoreInd.csv",sep="/"))$stkName
  colnames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars_CoreInd.csv",sep="/"))$stkName
}
if(AllExMH){
  rownames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars_AllExMH.csv",sep="/"))$stkName
  colnames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars_AllExMH.csv",sep="/"))$stkName
}
png(filename=paste(wcviCKDir, "/Figures/SpawnerCorrelation.png", sep=""), width=4, height=4.5, units="in", res=500)
par(xpd=TRUE)

#corrplot(corMat, method="circle", p.mat=corMat, insig="p-value", type="lower", mar=c(5,4,3,2))

# LW: I think if we add addCoef.col="black", tl.col = "black" to the corrplot function we should get the correlation coefficient numbers in all the cells. And tl.col="black" just makes the CU names black instead of red
#  LW: if we want to remove the cells with correlation values= 1 (e.g., a CU is perfectly correlated with itself), we can add diag=FALSE argument.

corrplot(corMat, method="circle", p.mat=corMat,  type="lower", addCoef.col="black", tl.col = "black", diag=FALSE, mar=c(0,0,0,0), tl.cex=0.9, cl.cex=0.7,  insig="p-value")
dev.off()
