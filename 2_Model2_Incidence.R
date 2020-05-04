#==========================================================================
#==========================================================================
# 2- Analysis of incidence (regitries data) and national estimation: Model2
#             Data= Lung cancer, female
#==========================================================================
#==========================================================================


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

load(paste0(repMain, "dataI.modified.RData"))


options(warn=-1)
rm(ky, kyc,bya, kage, kagec, distage, kamin, cage, cyear)
rm(dataI.ana, dataI.pred.fr, dataI.pred.zr, teI, respear, resdev,fit )
rm(sd.district, smooth.district, ranef.district,last.fix)
options(warn=1)

# Add centered AGE and YEAR 
#---------------------------
# Used in modelling to avoid potential numerical issues
cage=60
cyear=2000

dataI=dataI.modified
dataI$AGEC=dataI$AGE-cage
dataI$YEARC=dataI$YEAR-cyear

# Knots parametrisation 
ky=seq(1990,2010, by=5)
bya=10

# Data for analaysis analysis (district-level, districts covered by a registry)
dataI.ana = dataI[dataI$DISTRICT!="fr" & is.na(dataI$K)==F,]
dataI.ana$DISTRICT <- as.factor(dataI.ana$DISTRICT)

# Data for national predictions 
dataI.pred.fr = dataI[dataI$DISTRICT=="fr", ]
# Add cohort variable for visualization
dataI.pred.fr$COHORT <- dataI.pred.fr$YEAR-dataI.pred.fr$AGE

sink(paste0(repRes,"Model2_Incidence.lis"))
print(Sys.time())
print(sessionInfo())

# Data description 
#--------------------------------------
cat("\nDescription of data analayzed (number of cases)\n===================================\n")
print(tapply(dataI.ana$K, list(dataI.ana$YEAR, dataI.ana$DISTRICT), sum))

# Knots
#--------------------------------------
cat("\nKnots\n===================================\n")
#Position of interior knots from percentiles 0.01 up 90
distage <- aggregate(dataI.ana[ ,c("K" ,"PY")], by=list(AGE=dataI.ana$AGE), sum)
distage$AGE <- as.numeric(as.character(  distage$AGE))
distage$PCUM <-cumsum(distage$K)/sum(distage$K)
# Percentile 0.01, rounded at 10 (lower)
kamin <- distage[distage$PCUM >=0.01,][1,"AGE"]
kamin <- floor( kamin/bya)*bya
kage <- c(min(dataI.ana$AGE), seq(kamin,90,bya),99 )
kagec <- kage - cage

cat("\nAge: number and position of knots :\n")
cat( length(kage),", ", kage ,"\n\n") 
  
# Knots Year
kyc <- ky  - cyear
cat("\nYear: number and position of knots :\n")
cat( length(ky),", ", ky ,"\n\n") 
 

# Model 2
#--------------------------------------------
 
teI <- gam(K~te(AGEC, YEARC, k=c(length(kagec), length(kyc)),bs="cr")+ s(DISTRICT,bs="re") +  offset(log(PY)), knots=list(AGEC=kagec, YEARC=kyc), family=poisson, data=dataI.ana, method="REML")
 
cat("\nModel 2\n===================================\n")
cat("\nConvergence:\n")
print(teI$converged)
print(summary(teI))
cat("\nedf2:\n")
print(  sum(teI$edf2)) 
cat("\nSmoothing parameters\n")
print(teI$sp) 

# SD district random effect 
sd.district <- gam.vcomp(teI)
sd.district<-sd.district[match("s(DISTRICT)", row.names(sd.district)),1 ]
cat("\nSD district random effect:\n")
print(round(sd.district,3)) 


# Estimated district random effect 
# Note: district random effects are within the coefficient vector, after the MPS coefficients
# To identify the indices delimiting these coefficients, we use the smooth object 

smooth.district<- teI$smooth[[2]]
district.ranef <- exp(teI$coef[smooth.district$first.para:smooth.district$last.para])
names(district.ranef ) <- levels(dataI.ana$DISTRICT)
cat("\nEstimated random effects:\n")
print(round(district.ranef,2)) 

# Indice of the last fixed ceofficient (i.e. the last MPS coefficient) in the coef object (used below for prediction)
last.fix=smooth.district$first.para-1
  
 
# Prediction , national
#----------------------------------------------

dataI.pred.fr$K.pred <- rep(NA, nrow(dataI.pred.fr))

# Although DISTRICT not used for marginal prediction, it is required that
# the variable is in newdata, and with a known level 
# => temporarily set to "38" (one of the district)
dataI.pred.fr$DISTRICT<-as.factor(rep("38", nrow(dataI.pred.fr)))

dataI.pred.fr$K.pred <- exp(sd.district*sd.district/2)*predict.gam(teI, newdata= dataI.pred.fr, type="response", exclude="s(DISTRICT)")
dataI.pred.fr <- dataI.pred.fr[order(dataI.pred.fr$YEAR,dataI.pred.fr$AGE), ]
dataI.pred.fr$DISTRICT <- rep("fr",nrow(dataI.pred.fr))

 
 #Predicted incidence rate (la.rate) per 100000
 dataI.pred.fr$la.rate.pred <- 100000*dataI.pred.fr$K.pred/ dataI.pred.fr$PY
 
#Save national prediction by age and year
save(dataI.pred.fr, file=paste0(repRes,"dataI.pred.fr.RData"))
  
cat("\nPrediction, national level\n===================================\n")
cat("\nTotal number of incidence cases by year, predicted:\n")
print(round(tapply(dataI.pred.fr$K.pred, dataI.pred.fr$YEAR, sum),0))

cat("\nAverage annual number of cases 2011-2015:\n")
ind <- dataI.pred.fr$YEAR>=2011 & dataI.pred.fr$YEAR<=2015
print(round(sum(dataI.pred.fr$K.pred[ind]/5),0))


sink() 



