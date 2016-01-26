#AWKMT2 impact of choice of control group and some example null cases
library("plyr")      #, lib.loc="~/R/win-library/3.3")
library("dplyr", lib.loc="~/R/win-library/3.2")
library("readr", lib.loc="~/R/win-library/3.2")
library("rmarkdown", lib.loc="~/R/win-library/3.2")
library("survival", lib.loc="~/R/win-library/3.2")
library("controlTest", lib.loc="~/R/win-library/3.2")
library(bpcp);library(surv2sampleComp)
library(survRM2)                                                     #added Jan 10,2016
library(ggplot2);library(tidyr)
library("pander", lib.loc="~/R/win-library/3.2")
library("survTest2", lib.loc="~/R/win-library/3.2")              #Adaptively Weighted Kaplan-Meier tests
fp<-function(pvalues){            # pvalues<-a$by_strata[,1] or pvalues<-c(.6,.7)  
  df = 2*length(pvalues)
  1-pchisq( -2*sum(log(pvalues)), df, lower.tail=TRUE)
}


#-- sample data (pbc data in survival package)--#
D = pbc[1:312, c(2,3,4)] ; D$status=as.numeric(D$status==2) ; D$trt=as.numeric(D$trt==2)

tau=3000 ; 
nn=nrow(D)
st=sample(1:2, size=nn, replace=TRUE)
summary(D)
tau=max(D[D[,2]==1,1]); tau ;
indata=cbind(D, st)

a <- AWKMT2(indata, tau, stratified=TRUE)                               #Adaptive mean
AWKM_stats <-round(a$fisher,4)
V<-2                           #which statistic
print(paste(a$fisher[V],a$by_strata[1,V],a$by_strata[2,V],sep="  ") )

Su<-Surv(time=indata$time,event=indata$status)                     #KM rank test stratified 
MH<- survdiff(Su~indata$trt+strata(indata$st),rho=1)               
z<-sqrt(MH$chisq) 
p_Strat<-round( (1-pnorm(z)),4)                                    # one tail (1/2 of 2 tail)

Su<-Surv(time=indata$time,event=indata$status)                     #KM rank test stratified 
MH<- survdiff(Su~indata$trt+strata(indata$st),rho=1)               
z<-sqrt(MH$chisq) 
p_Strat<-round( (1-pnorm(z)),4)                                    # one tail (1/2 of 2 tail)




plot(survfit(Su~indata$trt+strata( indata$st) ),lty=c(1,2,1,2),col=c(1,2,3,4)) 
 
abline(v=tau,h=.5,col="lightgray")
title(main=paste( "Key: Wilcoxon 1/2 of 2 tail ", p_Strat," \n Adaptive Mean tau ",tau,"\n Fisher V1, V2   " ,
                  AWKM_stats[1] ,"   ", AWKM_stats[2] ))
 
# IMUC constructed data 
file="km_PFS_0_1.csv"
for (f in c(0,.1,.25,.333,.4,.5)){  
DD <- read_csv(file)
tau <- 3000
tau=max(DD[DD[,2]==1,1]); tau 
 
indata2 <- DD %>% transmute( time=OSTIME, status =DEATH, trt=as.numeric(TRT==1), st= as.factor(MGMT) )
# Dilute effect  by flipping group membership with probability f
nr<-nrow(indata2)

ran<-runif(nr)
for (III in 1: nr){
  if(ran[III] < f){
    ifelse(indata2$trt[III]==0,indata2$trt[III]<-1,indata2$trt[III]<-0)
  }
}

indata2<- as.data.frame(indata2)
a <- AWKMT2(indata2, tau, stratified=TRUE)                               #Adaptive mean
AWKM_stats <-round(a$fisher,4)      
print(paste(a$fisher,a$by_strata,sep="  ") )

Su<-Surv(time=indata2$time,event=indata2$status)                          #KM rank test stratified 
MH<- survdiff(Su~indata2$trt+strata(indata2$st),rho=1)               
z<-sqrt(MH$chisq) 
p_Strat<-round( (1-pnorm(z)),4)                                            # one tail
plot(survfit(Su~indata2$trt+strata( indata2$st) ),lty=c(1,2,1,2),col=c(1,2,3,4)) 
#Black and Green are strata 1 (solid)  Black  = Control  
#Red and Blue are strata 2   Red   = Control 3
 
abline(v=tau,h=.5,col="lightgray")
title(main=paste( "Key: (Black and Red are Control), (Strata 1 is solid, 2 is dashed)  ",
                  "\n Wilcoxon 1/2 of 2 tail ", p_Strat," \n Adaptive Mean tau ",tau,"\n Fisher V1, V2   " ,
                  AWKM_stats[1] , "      ",AWKM_stats[2] ),sub=paste(file, "dilute treatment with f =",f))

} 


