library(survival)
library(survRM2)
options(scipen=10)

nPerGroup <- 500

set.seed(1)

#Create Sample Trial Data
pid=seq(1, by=1, len=nPerGroup*2)
treatment=rep(1:2, nPerGroup)
survivaltime <- numeric(nPerGroup*2)
survivaltime[treatment==1] <- rweibull(nPerGroup, 1, scale0)
survivaltime[treatment==2] <- rweibull(nPerGroup, 1, scale1) 
sampletrial=data.frame(cbind(pid, treatment, survivaltime))
sampletrial$death<-ifelse(sampletrial$survivaltime <= 365, 1, 0)
sampletrial$deathtime<-ifelse(sampletrial$death == 1,sampletrial$survivaltime,365)
head(sampletrial, n=10)

#Create Indicator Variables For Alive/Dead at Final Analysis (365 Days)

sampletrial$deathdate=sampletrial$survivaltime
sampletrial$death[sampletrial$deathdate<365]=1
sampletrial$death[sampletrial$deathdate>=365]=0
head(sampletrial, n=10)

sampletrial$deathtime<-ifelse(sampletrial$death == 1,sampletrial$deathdate,365)
head(sampletrial, n=10)

table(sampletrial$death, sampletrial$treatment)

#Kaplan-Meier Curve

KM_final <- survfit(Surv(deathtime, death) ~ treatment, type="kaplan-meier", conf.type=c("none"), data=sampletrial)
plot(KM_final, main=expression(paste("Kaplan-Meier Estimate")), xlab="Time (Years)", ylab="Survival", lwd=2, col=1:2)

#Cox Model

coxmodel=coxph(Surv(deathtime, death) ~ treatment, data=sampletrial)
summary(coxmodel, conf.int=0.95)
hr=exp(summary(coxmodel)$coefficients[1])
hr
lcl=exp(summary(coxmodel)$coefficients[1]-1.96*summary(coxmodel)$coefficients[3])
lcl
ucl=exp(summary(coxmodel)$coefficients[1]+1.96*summary(coxmodel)$coefficients[3])
ucl
pvalue=summary(coxmodel)$coefficients[5]
pvalue

#RMST Analysis

time=sampletrial$deathtime
status=sampletrial$death
arm=sampletrial$treatment-1
rmst_results=rmst2(time, status, arm)
print(rmst_results)
plot(rmst_results)
rmst_diff=rmst_results$unadjusted.result[1]
rmst_diff
rmst_diff_lcl=rmst_results$unadjusted.result[4]
rmst_diff_lcl
rmst_diff_ucl=rmst_results$unadjusted.result[7]
rmst_diff_ucl
rmst_diff_pval=rmst_results$unadjusted.result[10]
rmst_diff_pval

#Store This Entire Code Sequence As Function

rm(trialresults)

trial_simulation <- function(nSims, nPerGroup, scale0, scale1, seed) {

set.seed(seed)

trial<-numeric(nSims)
hr <-numeric(nSims) 
hr_lcl <-numeric(nSims) 
hr_ucl <-numeric(nSims) 
hr_pvalue <-numeric(nSims) 
rmst_diff <-numeric(nSims) 
rmst_diff_lcl <-numeric(nSims) 
rmst_diff_ucl <-numeric(nSims) 
rmst_diff_pvalue <-numeric(nSims)
ni_rmst_1 <-numeric(nSims)  
ni_rmst_2 <-numeric(nSims)  
ni_rmst_3 <-numeric(nSims)  
ni_rmst_4 <-numeric(nSims)  

for(i in 1:nSims){

trial[i]=i

#Create Sample Trial Data
pid=seq(1, by=1, len=nPerGroup*2)
treatment=rep(1:2, nPerGroup)
survivaltime <- numeric(nPerGroup*2)
survivaltime[treatment==1] <- rweibull(nPerGroup, 1, scale0)
survivaltime[treatment==2] <- rweibull(nPerGroup, 1, scale1) 
sampletrial=data.frame(cbind(pid, treatment, survivaltime))
sampletrial$death<-ifelse(sampletrial$survivaltime <= 365, 1, 0)
sampletrial$deathtime<-ifelse(sampletrial$death == 1,sampletrial$survivaltime,365)

#Cox Model
coxmodel=coxph(Surv(deathtime, death) ~ treatment, data=sampletrial)
hr[i]=exp(summary(coxmodel)$coefficients[1])
hr_lcl[i]=exp(summary(coxmodel)$coefficients[1]-1.96*summary(coxmodel)$coefficients[3])
hr_ucl[i]=exp(summary(coxmodel)$coefficients[1]+1.96*summary(coxmodel)$coefficients[3])
hr_pvalue[i]=summary(coxmodel)$coefficients[5]

#RMST Analysis
time=sampletrial$deathtime
status=sampletrial$death
arm=sampletrial$treatment-1
rmst_results=rmst2(time, status, arm)
rmst_diff[i]=rmst_results$unadjusted.result[1]
rmst_diff_lcl[i]=rmst_results$unadjusted.result[4]
rmst_diff_ucl[i]=rmst_results$unadjusted.result[7]
rmst_diff_pvalue[i]=rmst_results$unadjusted.result[10]
ni_rmst_1[i]<-ifelse(rmst_diff_lcl[i] > - 7 , 1, 0)
ni_rmst_2[i]<-ifelse(rmst_diff_lcl[i] > - 14 , 1, 0)
ni_rmst_3[i]<-ifelse(rmst_diff_lcl[i] > - 21 , 1, 0)
ni_rmst_4[i]<-ifelse(rmst_diff_lcl[i] > - 28 , 1, 0)

}

trialresults=as.data.frame(cbind(trial, 
hr, hr_lcl, hr_ucl, hr_pvalue, 
rmst_diff, rmst_diff_lcl, rmst_diff_ucl, rmst_diff_pvalue,
ni_rmst_1, ni_rmst_2, ni_rmst_3, ni_rmst_4))

return(trialresults)

}

trialresults1=trial_simulation(nSims=1000, nPerGroup=500, scale0=10000, scale1=10000, seed=1)
table(trialresults1$ni_rmst_1)
table(trialresults1$ni_rmst_2)
table(trialresults1$ni_rmst_3)
table(trialresults1$ni_rmst_4)

trialresults2=trial_simulation(nSims=1000, nPerGroup=500, scale0=5000, scale1=5000, seed=1)
table(trialresults2$ni_rmst_1)
table(trialresults2$ni_rmst_2)
table(trialresults2$ni_rmst_3)
table(trialresults2$ni_rmst_4)

trialresults3=trial_simulation(nSims=1000, nPerGroup=500, scale0=3000, scale1=3000, seed=1)
table(trialresults3$ni_rmst_1)
table(trialresults3$ni_rmst_2)
table(trialresults3$ni_rmst_3)
table(trialresults3$ni_rmst_4)

trialresults4=trial_simulation(nSims=1000, nPerGroup=500, scale0=2000, scale1=2000, seed=1)
table(trialresults4$ni_rmst_1)
table(trialresults4$ni_rmst_2)
table(trialresults4$ni_rmst_3)
table(trialresults4$ni_rmst_4)

trialresults5=trial_simulation(nSims=1000, nPerGroup=500, scale0=1000, scale1=1000, seed=1)
table(trialresults5$ni_rmst_1)
table(trialresults5$ni_rmst_2)
table(trialresults5$ni_rmst_3)
table(trialresults5$ni_rmst_4)








