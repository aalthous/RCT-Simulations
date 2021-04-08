library(survival)
library(truncnorm)
library(BuyseTest)
options(scipen=10)

trial_simulation <- function(nSims, nPerGroup, scale1a, scale2a, scale1b, scale2b, seed) {

set.seed(seed)

trial<-numeric(nSims)
hr<-numeric(nSims) 
hr_lcl<-numeric(nSims) 
hr_ucl<-numeric(nSims) 
hr_pvalue<-numeric(nSims) 
hr_success<-numeric(nSims)
wr<-numeric(nSims)
wr_lcl<-numeric(nSims)
wr_ucl<-numeric(nSims)
wr_pvalue<-numeric(nSims)
wr_success<-numeric(nSims)

for(i in 1:nSims){

trial[i]=i

#Create Sample Trial Data
pid=seq(1, by=1, len=nPerGroup*2)
risk=rtruncnorm(nPerGroup, a=0, b=1, mean=0.5, sd=0.1)
hosptimes1 <- (1-risk)*rweibull(nPerGroup,1,scale1a)
hosptimes2 <- (1-risk)*rweibull(nPerGroup,1,scale2a)
survtimes1 <- (1-risk)*rweibull(nPerGroup,1,scale1b)
survtimes2 <- (1-risk)*rweibull(nPerGroup,1,scale2b)
treatment=rep(0:1, nPerGroup)
hosptime <- numeric(nPerGroup*2)
hosptime <- ifelse(treatment==0, hosptimes1, hosptimes2)
survivaltime <- numeric(nPerGroup*2)
survivaltime<- ifelse(treatment==0, survtimes1, survtimes2)

death <- ifelse(survivaltime <= 365, 1, 0)
deathtime<-ifelse(death == 1, survivaltime, 365)

hfhosp<-ifelse(hosptime <= 365 & hosptime <=survivaltime, 1, 0)
hfhosptime<-ifelse(hfhosp == 1, hosptime, 365)

hfhosp<-ifelse(death == 1 & hosptime>survivaltime, 0, hfhosp)
hfhosptime<-ifelse(death==1 & hosptime>survivaltime, deathtime, hfhosptime)

death_hfh<-ifelse(death==1, 1, hfhosp)
death_hfh_time<-ifelse(death_hfh==1, hfhosptime, deathtime)

sampletrial=data.frame(cbind(pid, treatment, hfhosp, hfhosptime, death, deathtime, death_hfh, death_hfh_time))

#Cox Model for Death & HFH
coxmodel=coxph(Surv(death_hfh_time, death_hfh) ~ treatment, data=sampletrial)
hr[i]=exp(summary(coxmodel)$coefficients[1])
hr_lcl[i]=exp(summary(coxmodel)$coefficients[1]-1.96*summary(coxmodel)$coefficients[3])
hr_ucl[i]=exp(summary(coxmodel)$coefficients[1]+1.96*summary(coxmodel)$coefficients[3])
hr_pvalue[i]=summary(coxmodel)$coefficients[5]
hr_success[i]<-ifelse(hr[i]<1 & hr_pvalue[i] < 0.05, 1, 0)

# Win Ratio for Death & HFH
Setup <- treatment ~ tte(deathtime, status="death") + tte(hfhosptime, status="hfhosp")
BT <- BuyseTest(Setup, data=sampletrial, trace=0)
winratio<-summary(BT, statistic="WinRatio", print=FALSE)
wr[i]<-winratio$table.print$Delta[2]
wr_lcl[i]<-winratio$table$CIinf.Delta[3]
wr_ucl[i]<-winratio$table$CIsup.Delta[3]
wr_pvalue[i]<-winratio$table$p.value[3]
wr_success[i]<-ifelse(wr[i]>1 & wr_pvalue[i] < 0.05, 1, 0)
}

trialresults=as.data.frame(cbind(trial, 
hr, hr_lcl, hr_ucl, hr_pvalue, hr_success,
wr, wr_lcl, wr_ucl, wr_pvalue, wr_success))

return(trialresults)

}

trialresults1=trial_simulation(nSims=100, nPerGroup=500, scale1a=750, scale2a=1000, scale1b=1500, scale2b=2000, seed=1)
head(trialresults1,n=10)
table(trialresults1$hr_success)
table(trialresults1$wr_success)
