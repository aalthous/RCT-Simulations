set.seed(1)
library(survival) 
nPerGroup <- 500 
recruitPerDay <- 1
weibull_shape1 <- 0.5 
weibull_shape2 <- 0.5 
weibull_scale1 <- 1200 
weibull_scale2 <- 2000 
interim_analysis_dates <- c(720, 1080, 1440)
final_analysis_date <- 1800 

#Create One Sample Trial Dataset
pid=seq(1, by=1, len=nPerGroup*2) 
treatment=rep(1:2, nPerGroup)
recruitdate=floor((pid+1)/recruitPerDay) 
survivaltime <- numeric(nPerGroup*2) 
survivaltime[treatment==1] <- rweibull(nPerGroup, weibull_shape1, weibull_scale1)
survivaltime[treatment==2] <- rweibull(nPerGroup, weibull_shape2, weibull_scale2) 

deathdate=recruitdate+survivaltime

death1 <- numeric(nPerGroup*2)
death1[deathdate<interim_analysis_dates[1]]=1
death1[deathdate>=interim_analysis_dates[1]]=0

deathtime1<- 
ifelse(death1 == 1,
deathdate-recruitdate,
interim_analysis_dates[1]-recruitdate)

death2 <- numeric(nPerGroup*2)
death2[deathdate<interim_analysis_dates[2]]=1
death2[deathdate>=interim_analysis_dates[2]]=0

deathtime2<- 
ifelse(death2 == 1,
deathdate-recruitdate,
interim_analysis_dates[2]-recruitdate)

death3 <- numeric(nPerGroup*2)
death3[deathdate<interim_analysis_dates[3]]=1
death3[deathdate>=interim_analysis_dates[3]]=0

deathtime3<- 
ifelse(death3 == 1,
deathdate-recruitdate,
interim_analysis_dates[3]-recruitdate)

death_final <- numeric(nPerGroup*2)
death_final[deathdate<final_analysis_date]=1
death_final[deathdate>=final_analysis_date]=0

deathtime_final<- 
ifelse(death == 1,
deathdate-recruitdate,
final_analysis_date-recruitdate) 

sampletrial <- data.frame(cbind(pid, treatment, recruitdate, death1, deathtime1, death2, deathtime2, death3, deathtime3, death_final, deathtime_final))

#Kaplan-Meier Curve for Sample Trial Data

KM_final <- survfit(Surv(deathtime_final, death_final) ~ treatment, type="kaplan-meier", conf.type=c("none"), data=sampletrial)
plot(KM_final, main=expression(paste("Kaplan-Meier Estimate")), xlab="Time (Years)", ylab="Survival", lwd=2, col=1:2)

coxmodel=coxph(Surv(deathtime, death) ~ treatment, data=sampletrial)
summary(coxmodel)

nLooks<-4
analyses_scheduled<-(c(0.40, 0.60, 0.80, 1))
efficacy_thresholds<-numeric(4)
design <- getDesignGroupSequential(sided=1, alpha=0.05, informationRates=analyses_scheduled, typeOfDesign = "asOF")
for(j in 1:nLooks){
efficacy_thresholds[j] = design$stageLevels[j]
}
efficacy_thresholds

# Simulation Parameters
nSims <- 1000
trialnum <- numeric(nSims)
hr <- data.frame(matrix(ncol = nLooks, nrow = nSims))
lcl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
ucl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
pvalue <- data.frame(matrix(ncol = nLooks, nrow = nSims))
success <- data.frame(matrix(ncol = nLooks, nrow = nSims))
overall_success <-numeric(nSims)

#provide column names
colnames(hr) <- sprintf("hr_%d", (1:nLooks))
colnames(lcl) <- sprintf("lcl_%d", (1:nLooks))
colnames(ucl) <- sprintf("ucl_%d", (1:nLooks))
colnames(pvalue) <- sprintf("pvalue_%d", (1:nLooks))
colnames(success) <- sprintf("success_%d", (1:nLooks))
set.seed(1) # this sets the random seed for your results to be reproducible

for(i in 1:nSims){

trialnum[i]=i 

pid=seq(1, by=1, len=nPerGroup*2) 
treatment=rep(1:2, nPerGroup)
recruitdate=floor((pid+1)/recruitPerDay) 
survivaltime <- numeric(nPerGroup*2) 
survivaltime[treatment==1] <- rweibull(nPerGroup, weibull_shape1, weibull_scale1)
survivaltime[treatment==2] <- rweibull(nPerGroup, weibull_shape2, weibull_scale2) 
deathdate=recruitdate+survivaltime

death1 <- numeric(nPerGroup*2)
death1[deathdate<interim_analysis_dates[1]]=1
death1[deathdate>=interim_analysis_dates[1]]=0

deathtime1<- 
ifelse(death1 == 1,
deathdate-recruitdate,
interim_analysis_dates[1]-recruitdate)

death2 <- numeric(nPerGroup*2)
death2[deathdate<interim_analysis_dates[2]]=1
death2[deathdate>=interim_analysis_dates[2]]=0

deathtime2<- 
ifelse(death2 == 1,
deathdate-recruitdate,
interim_analysis_dates[2]-recruitdate)

death3 <- numeric(nPerGroup*2)
death3[deathdate<interim_analysis_dates[3]]=1
death3[deathdate>=interim_analysis_dates[3]]=0

deathtime3<- 
ifelse(death3 == 1,
deathdate-recruitdate,
interim_analysis_dates[3]-recruitdate)

death_final <- numeric(nPerGroup*2)
death_final[deathdate<final_analysis_date]=1
death_final[deathdate>=final_analysis_date]=0

deathtime_final<- 
ifelse(death == 1,
deathdate-recruitdate,
final_analysis_date-recruitdate) 

trialdata <- data.frame(cbind(pid, treatment, recruitdate,
death1, deathtime1, 
death2, deathtime2, 
death3, deathtime3, 
death_final, deathtime_final))

analysisdata1 <- subset(trialdata, recruitdate<interim_analysis_dates[1])
coxmodel1=coxph(Surv(deathtime1, death1) ~ treatment, data=analysisdata1)
hr[i,1]=exp(summary(coxmodel1)$coefficients[1])
lcl[i,1]=exp(summary(coxmodel1)$coefficients[1]-1.96*summary(coxmodel1)$coefficients[3]) 
ucl[i,1]=exp(summary(coxmodel1)$coefficients[1]+1.96*summary(coxmodel1)$coefficients[3]) 
pvalue[i,1]=summary(coxmodel1)$coefficients[5] 
success[i,1]=ifelse(hr[i,1]<1 & pvalue[i,1]<efficacy_thresholds[1], 1, 0) 

analysisdata2 <- subset(trialdata, recruitdate<interim_analysis_dates[2])
coxmodel2=coxph(Surv(deathtime2, death2) ~ treatment, data=analysisdata2)
hr[i,2]=exp(summary(coxmodel2)$coefficients[1])
lcl[i,2]=exp(summary(coxmodel2)$coefficients[1]-1.96*summary(coxmodel2)$coefficients[3]) 
ucl[i,2]=exp(summary(coxmodel2)$coefficients[1]+1.96*summary(coxmodel2)$coefficients[3]) 
pvalue[i,2]=summary(coxmodel2)$coefficients[5] 
success[i,2]=ifelse(hr[i,2]<1 & pvalue[i,2]<efficacy_thresholds[2], 1, 0) 

analysisdata3 <- subset(trialdata, recruitdate<interim_analysis_dates[3])
coxmodel3=coxph(Surv(deathtime3, death3) ~ treatment, data=analysisdata3)
hr[i,3]=exp(summary(coxmodel3)$coefficients[1])
lcl[i,3]=exp(summary(coxmodel3)$coefficients[1]-1.96*summary(coxmodel3)$coefficients[3]) 
ucl[i,3]=exp(summary(coxmodel3)$coefficients[1]+1.96*summary(coxmodel3)$coefficients[3]) 
pvalue[i,3]=summary(coxmodel3)$coefficients[5] 
success[i,3]=ifelse(hr[i,3]<1 & pvalue[i,3]<efficacy_thresholds[3], 1, 0) 

coxmodel_final=coxph(Surv(deathtime3, death3) ~ treatment, data=trialdata)
hr[i,4]=exp(summary(coxmodel_final)$coefficients[1])
lcl[i,4]=exp(summary(coxmodel_final)$coefficients[1]-1.96*summary(coxmodel_final)$coefficients[3]) 
ucl[i,4]=exp(summary(coxmodel_final)$coefficients[1]+1.96*summary(coxmodel_final)$coefficients[3]) 
pvalue[i,4]=summary(coxmodel_final)$coefficients[5] 
success[i,4]=ifelse(hr[i,4]<1 & pvalue[i,4]<efficacy_thresholds[4], 1, 0) 

overall_success[i] <- 0

for (j in 1:nLooks)
{
if(success[i,j]==1)
{
overall_success[i] <- 1 
}
}

}

simulation_results <- data.frame(cbind(trialnum,
hr[1], lcl[1], ucl[1], pvalue[1], success[1],
hr[2], lcl[2], ucl[2], pvalue[2], success[2],
hr[3], lcl[3], ucl[3], pvalue[3], success[3],
hr[4], lcl[4], ucl[4], pvalue[4], success[4],
overall_success))
head(simulation_results, n=10)

table(overall_success)
table(simulation_results$success_1, overall_success)
table(simulation_results$success_2, overall_success)
table(simulation_results$success_3, overall_success)
table(simulation_results$success_4, overall_success)
