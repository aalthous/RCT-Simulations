library(truncnorm)
library(survival)
library(BuyseTest)
options(scipen=10)

trial_simulation <- function(nSims, nPerGroup, maxtime, scale0, scale1, trt_lvef, seed) {

set.seed(seed)

trial<-numeric(nSims)
wr<-numeric(nSims)
wr_lcl<-numeric(nSims)
wr_ucl<-numeric(nSims)
wr_pvalue<-numeric(nSims)
wr_success<-numeric(nSims)

for(i in 1:nSims){

trial[i]=i

pid <- seq(1, by=1, len=nPerGroup*2)

treatment=rep(0:1, nPerGroup)

bl_lvef <- rtruncnorm(nPerGroup*2, a=0, b=0.35, mean=0.25, sd=0.1)

survtimes0 <- rweibull(nPerGroup*2,1,scale0)
survtimes1 <- rweibull(nPerGroup*2,1,scale1)

survivaltime <- numeric(nPerGroup*2)
survivaltime <- ifelse(treatment==0, survtimes0, survtimes1)

death <- ifelse(survivaltime <= maxtime, 1, 0)
deathtime<-ifelse(death == 1, survivaltime, maxtime)

ch_lvef <- rtruncnorm(nPerGroup*2, a=0, b=0.30, mean=0.15, sd=0.1) + treatment*trt_lvef

fu_lvef <- ifelse(death == 1, NA, bl_lvef + ch_lvef)

sampletrial <- data.frame(cbind(pid, treatment, bl_lvef, death, deathtime, fu_lvef))

# Win Ratio for Death & LVEF
Setup <- treatment ~ tte(deathtime, status="death") + cont(fu_lvef)
BT <- BuyseTest(Setup, data=sampletrial, trace=0)
winratio <- summary(BT, statistic="WinRatio", print=FALSE)
wr[i] <- winratio$table.print$Delta[2]
wr_lcl[i] <- winratio$table$CIinf.Delta[3]
wr_ucl[i] <- winratio$table$CIsup.Delta[3]
wr_pvalue[i] <- winratio$table$p.value[3]
wr_success[i] <- ifelse(wr[i] > 1 & wr_pvalue[i] < 0.05, 1, 0)
}

trialresults <- as.data.frame(cbind(trial,
wr, wr_lcl, wr_ucl, wr_pvalue, wr_success))

return(trialresults)

}

#maxtime=180 -> assuming LVEF assessed at 6 months
#scale0=2000 -> assigns a 6-month event rate ~8.6% for patients with treatment=0 based on pweibull(180,1,2000)
#scale1=5000 -> assigns a 6-month event rate ~3.5% for patients with treatment=1 based on pweibull(180,1,5000)
#trt_lvef=0.05 -> assigns the "treatment effect" of LVEF improvement being +0.05 for treatment=1 vs treatment=0

trialresults1 <- trial_simulation(nSims=1000, nPerGroup=90, maxtime=180, scale0=2000, scale1=5000, trt_lvef=0.05, seed=1)
head(trialresults1,n=10)
table(trialresults1$wr_success)
