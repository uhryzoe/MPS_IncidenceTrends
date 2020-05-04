#==========================================================================
#==========================================================================
# 3- Example of visualization : predicted mortality 
#             Data analyzed= Lung cancer female
#==========================================================================
#==========================================================================


#============================================
# Directories (main: name to be changed) 
# NOTE: creates a /res subdirectory in the main directory
#======================================================


# ***!! Name to changed by user !!***
repMain <- "D/Temp/"
setwd(repMain)

repRes <- "./res/" 
# Create a /res  subdirectoty
dir.create(repRes)

options(warn=-1)
rm(dataM.pred)
rm(lage,lcohortccol, cpch)
options(warn=1)

load(paste0(repRes,"dataM.pred.fr.Rdata"))


# 3D-plot vizulation of prediction (similar Figure 1)
#=======================================================================
res <- dataM.pred.fr
res <- res[order(res$YEAR, res$AGE), ]

pdf(file=paste(repRes, "/Mortality-Plot3D.pdf", sep=""), paper="special" ,width=12/2.54, height=10/2.54)
par(oma=c(0.01,0.01,0.01,0.01), mai=c(0.8,0.8,0.5,0.2))
ind <- res$YEAR%in%seq(1976,2018, by=2) &res$AGE%in%seq(20,96,2) 
b<-  reshape(res[ind, c("mu.rate.pred","YEAR","AGE")], timevar = "YEAR", idvar="AGE", direction = "wide")
b <- b[ , -1]
persp(x=seq(20,96,2), y=seq(1976,2018, by=2),z=as.matrix(b),ylim=c(1975,2018), ticktype="detailed",xlab="\nAge",ylab="\n\nYear",zlab="\n\nMortality rate",theta=30,phi=20,expand=0.9,d=3,lwd=0.5,cex.axis=0.8, cex.lab=0.8, col="lightblue")
mtext("Female lung cancer mortality", cex=0.8)
dev.off()


# Trends in mortality rate by age, yyear or cohort (similar to the second column of figure 3)
#=======================================================================
res <- dataM.pred.fr
res <- res[order(res$YEAR, res$AGE), ]

lage= seq(40,80, by=10)
lcohort=seq(1910,1970, by=10)
lyear=c(seq(1990,2010, by=10),2015,2018)
ccol=c(1:2,"chartreuse4",4,"darkcyan",3,6 )
cpch=c(1:2,8,4:5,3,6)

pdf(file=paste(repRes, "/Mortality-Visualization.pdf", sep=""),  paper="special" ,width=20/2.54, height=17/2.54)
par(mfrow=c(2,2), oma=c(0.01,0.01,2,0.01), mai=c(0.8,0.6,0.4,0.4), mgp=c(2.5,1,0))

# a) Plot transversal age-curve (for different years) 
#-----------------------------------------------------

a <- res[res$YEAR%in%lyear, ]
ymax <- 1.2*max(a$mu.rate.pred,na.rm=T)
ymin <-0.6*min(a$mu.rate.pred,na.rm=T)
i=0
for (year in lyear)
{
i<- i+1
xa <- a[a$YEAR==year & a$AGE>=30 & a $AGE<=90 , ]
indp <- xa$AGE%in%seq(40,80, by=10)
if (i==1) { plot(xa$AGE, xa$mu.rate.pred, type="l",xlim=c(30,90),ylim=c(ymin,ymax),  xlab="Age",ylab="", cex.axis=1.2, cex.lab=1.2) }
if (year<2018){lines(xa$AGE, xa$mu.rate.pred,col=ccol[i], lty=1)}
if (year==2018){lines(xa$AGE, xa$mu.rate.pred,col=ccol[i], lty=2)}
points(xa$AGE[indp], xa$mu.rate.pred[indp],col=ccol[i], pch=cpch[i])
}
legend("topleft", legend=as.character(lyear), col=ccol, pch=cpch, bty="n", horiz=F, cex=1.1, title="Year")
mtext("Age-specific rate for different years, Mortality (a)", side=3, line=0.7, cex=0.8)
mtext("Mortality rate", side=2, line=2, cex=0.9)



#  b) Plot longitudinal age-curve (for different cohorts)
#-----------------------------------------------------

a <- res[res$COHORT%in%lcohort, ]
ymax <- 1.2*max(a$mu.rate.pred,na.rm=T)
ymin <-0.6*min(a$mu.rate.pred,na.rm=T)
i=0
for (cohort in lcohort)
{
i<- i+1
xa <- a[a$COHORT==cohort & a$AGE>=30 & a $AGE<=90, ]
ind <- (xa$YEAR <=2015)
indp <- xa$AGE%in%seq(20,90, by=10)
if (i==1) { plot(xa$AGE[ind], xa$mu.rate.pred[ind], type="l",xlim=c(30,90),ylim=c(ymin,ymax), , xlab="Age",ylab="", cex.axis=1.2, cex.lab=1.2) }
lines(xa$AGE, xa$mu.rate.pred,col=ccol[i], lty=2)
lines(xa$AGE[ind], xa$mu.rate.pred[ind],col=ccol[i])
points(xa$AGE[indp], xa$mu.rate.pred[indp],col=ccol[i], pch=cpch[i])
}
legend("topleft",  legend=as.character(lcohort), col=ccol, pch=cpch, bty="n", horiz=F, cex=1.1, title="Cohort")
mtext("Age-specific rate for different cohorts, Mortality (b)", side=3, line=0.7, cex=0.8)
mtext("Mortality rate", side=2, line=2, cex=0.9)


#  c) Plot trend by year for different age 
#-----------------------------------------------------
a <- res[res$AGE%in%lage & res$YEAR>=1990, ]
ymax <- 1.2*max(a$mu.rate.pred,na.rm=T)
ymin <-0.6*min(a$mu.rate.pred,na.rm=T)
i=0
for (age in lage)
{
i<- i+1
xa <- a[a$AGE==age, ]
ind <- (xa$YEAR <=2015)
indp <- xa$YEAR%in%seq(1980,2010,5)
if (i==1) { plot(xa$YEAR[ind], xa$mu.rate.pred[ind], type="l",xlim=c(1990,2018),ylim=c(ymin,ymax), log="y", xlab="Year",ylab="", cex.axis=1.2, cex.lab=1.2) }
lines(xa$YEAR, xa$mu.rate.pred,col=ccol[i], lty=2)
lines(xa$YEAR[ind], xa$mu.rate.pred[ind],col=ccol[i])
points(xa$YEAR[indp], xa$mu.rate.pred[indp],col=ccol[i], pch=cpch[i])
}
legend("bottomleft", c("Age:", lage), col=c(-1,ccol), pch=c(-1,cpch), bty="n", horiz=T, cex=1.1)
mtext("Trends by year for different ages (log-scale), Mortality (c)", side=3, line=0.7, cex=0.8)
mtext("Mortality rate", side=2, line=2, cex=0.9)


#  d) Plot trend by cohort for different age 
#-----------------------------------------------------
a <- res[res$AGE%in%lage & res$COHORT>=1910, ]
ymax <- 1.2*max(a$mu.rate.pred,na.rm=T)
ymin <-0.6*min(a$mu.rate.pred,na.rm=T)
i=0
for (age in lage)
{
i<- i+1
xa <- a[a$AGE==age, ]
ind <- (xa$YEAR <=2015)
indp <- xa$COHORT%in%lcohort
if (i==1) { plot(xa$COHORT[ind], xa$mu.rate.pred[ind], type="l",xlim=c(1910,1980),ylim=c(ymin,ymax), log="y", xlab="Cohort",ylab="", cex.axis=1.2, cex.lab=1.2) }
lines(xa$COHORT, xa$mu.rate.pred,col=ccol[i], lty=2)
lines(xa$COHORT[ind], xa$mu.rate.pred[ind],col=ccol[i])
points(xa$COHORT[indp], xa$mu.rate.pred[indp],col=ccol[i], pch=cpch[i])
}
legend("bottomleft", c("Age:", lage), col=c(-1,ccol), pch=c(-1,cpch), bty="n", horiz=T, cex=1.1)
mtext("Trends by cohort for different ages (log-scale), Mortality (d)", side=3, line=0.7, cex=0.8)
mtext("Mortality rate", side=2, line=2, cex=0.9)

mtext("Female lung cancer mortality", outer=T, line=0.1)

dev.off()


