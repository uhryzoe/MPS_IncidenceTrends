#===========================================================
#===========================================================
# 1- Analysis of national mortality data: Model 1
#             Data= Lung cancer, female
#===========================================================
#===========================================================

# Libraries
#============================================
library(mgcv)

#============================================
# Directories (main: name to be changed) and data load
# NOTE: creates a /res subdirectory in the main directory
#======================================================

# ***!! Name to changed by user !!***
repMain <- "D/Temp/"
setwd(repMain)

repRes <- "./res/" 
# Create a /res  subdirectoty
dir.create(repRes)

load(paste0(repMain, "dataM.RData"))


options(warn=-1)
rm(ky, kyc,bya, kage, kagec, distage, kamin, cage, cyear)
rm(dataM.ana, dataM.pred.fr, teM, respear, resdev,fit )
options(warn=1)

# Add centered AGE and YEAR 
#---------------------------
# Used in modelling to avoid potential numerical issues
cage=60
cyear=2000

dataM$AGEC=dataM$AGE-cage
dataM$YEARC=dataM$YEAR-cyear

# Knots parametrisation 
ky=seq(1980,2010, by=5)
bya=10

# Data for analysis (years 2016-2018 excluded)
dataM.ana = dataM[is.na(dataM$DEATH)==F,]

# Data for prediction (years 2016-2018 included)
dataM.pred.fr = dataM
# Add cohort variable  for visualization
dataM.pred.fr$COHORT <- dataM.pred.fr$YEAR-dataM.pred.fr$AGE

sink(paste0(repRes,"Model1_Mortality.lis"))
print(Sys.time())
print(sessionInfo())

# Data description 
#--------------------------------------
cat("\nDescription of data analyzed (number of death)\n===================================\n")
print(tapply(dataM.ana$DEATH, dataM.ana$YEAR, sum))


# Knots
#--------------------------------------
cat("\nKnots\n===================================\n")
 #Position of interior knots from percentiles 0.01 up 90
 distage <- aggregate(dataM.ana[ ,c("DEATH" ,"PY")], by=list(AGE=dataM.ana$AGE), sum)
 distage$AGE <- as.numeric(as.character(  distage$AGE))
 distage$PCUM <-cumsum(distage$DEATH)/sum(distage$DEATH)
 # Percentile 0.01, rounded at 10 (lower)
 kamin <- distage[distage$PCUM >=0.01,][1,"AGE"]
 kamin <- floor( kamin/bya)*bya
 kage <- c(min(dataM.ana$AGE), seq(kamin,90,bya),99 )
 kagec <- kage - cage

 cat("\nAge: number and position of knots :\n")
 cat( length(kage),", ", kage ,"\n\n") 
   
 # Knots Year
 kyc <- ky  - cyear
 cat("\nYear: number and position of knots :\n")
 cat( length(ky),", ", ky ,"\n\n") 
 

# Model 1
#--------------------------------------------
 
 teM <- gam(DEATH~ te(AGEC,YEARC,k=c(length(kagec),length(kyc)),bs="cr")+ offset(log(PY)),   knots=list(AGEC=kagec, YEARC=kyc), family=poisson, data=dataM.ana, method="REML")
 
 cat("\nModel 1\n===================================\n")
 cat("\nConvergence:\n")
 print(teM$converged)
 print(summary(teM))
 cat("\nedf2:\n")
 print(  sum(teM$edf2)) 
 cat("\nSmoothing parameters\n")
 print(teM$sp) 

  
# Prediction 
#----------------------------------------------
 
 cat("\nPrediction\n===================================\n")
 dataM.pred.fr$DEATH.pred<- predict.gam(teM, dataM.pred.fr, type = "response")
 
 #Observed and predicted mortality rate (mu.rate) per 100000
 dataM.pred.fr$mu.rate.obs <- 100000*dataM.pred.fr$DEATH/ dataM.pred.fr$PY
 dataM.pred.fr$mu.rate.pred <- 100000*dataM.pred.fr$DEATH.pred/ dataM.pred.fr$PY
 
 dataM.pred.fr<- dataM.pred.fr[order(dataM.pred.fr$SEXE,dataM.pred.fr$YEAR,dataM.pred.fr$AGE), ]

 #Save national prediction by age and year
 save(dataM.pred.fr, file=paste0(repRes,"dataM.pred.fr.RData"))
 
 cat("\nTotal number of death by year, observed:\n")
 print(round(tapply(dataM.pred.fr$DEATH, dataM.pred.fr$YEAR, sum),0))
 cat("\nTotal number of death by year, predicted:\n")
 print(round(tapply(dataM.pred.fr$DEATH.pred, dataM.pred.fr$YEAR, sum),0))
 
 cat("\nRelative error on the total number of death\n(predicted-observed)/observed x 100:\n")
 pred.tot<- tapply(dataM.pred.fr$DEATH.pred, dataM.pred.fr$YEAR, sum)
 obs.tot <-tapply(dataM.pred.fr$DEATH, dataM.pred.fr$YEAR, sum)
 print(round(100*((pred.tot-obs.tot)/obs.tot),1))
 
 sink()
 
 