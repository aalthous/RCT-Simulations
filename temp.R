library(MASS)
library(lme4)
library(lmerTest)
options(scipen=10)

TrialSimulation <- function(seed, nSims, nSites, nPerSite, nTimePts, site_var, muCtrl, muIntv, stdevs, mvncorr, weibull_param)
{

set.seed(seed)

trial <- numeric(nSims)
sd_site <- numeric(nSims)
ICC <- numeric(nSims)
b_1 <- numeric(nSims) 
lcl_1 <- numeric(nSims) 
ucl_1 <- numeric(nSims) 
pval_1 <- numeric(nSims) 
succ_1 <- numeric(nSims)
b_2 <- numeric(nSims) 
lcl_2 <- numeric(nSims) 
ucl_2 <- numeric(nSims) 
pval_2 <- numeric(nSims) 
succ_2 <- numeric(nSims)
b_3 <- numeric(nSims) 
lcl_3 <- numeric(nSims) 
ucl_3 <- numeric(nSims) 
pval_3 <- numeric(nSims) 
succ_3 <- numeric(nSims)

nTotal <- nSites*nPerSite

for(i in 1:nSims){
  
trial[i]=i

Site<-seq(1, by=1, len=nSites)
Site_Assign<-rep(0:1, nSites/2)
Site_Assignments<-as.data.frame(cbind(Site, Site_Assign))

Site_Effect<-rnorm(nSites, 0, site_var)
Site_Effects<-as.data.frame(cbind(Site, Site_Effect))

Site<-floor(runif(nTotal, 1, nSites+1))
z=table(Site)
sd_site[i]=sd(z)

Subject <- factor(1:(nTotal))
SurvivalTime <- rweibull(nTotal, 1, weibull_param)
Subjects <- as.data.frame(cbind(Subject, Site, SurvivalTime))

Setup1 = merge(x=Site_Assignments, y=Subjects, by="Site", all=TRUE)
Setup2 = merge(x=Setup1, y=Site_Effects, by="Site", all=TRUE)
Setup3 = Setup2[order(Setup2$Subject),]

Time <- factor(1:nTimePts, labels = c("PRE", "Time1", "Time2", "Time3"))
Setup4 <- expand.grid(Time, Subject)
names(Setup4) <- c("Time", "Subject")

Time2 <- rep(0:3, nTotal)
Setup5 <- as.data.frame(cbind(Setup4, Time2))

Setup6 = merge(x=Setup3, y=Setup5, by="Subject", all=TRUE)

stdiff <- mean(stdevs)
ones <- rep(1, nTimePts)
A <- stdevs^2 %o% ones
B <- (A + t(A) + (stdiff^2)*(diag(nTimePts) - ones %o% ones))*mvncorr

collaboRATE_ctrl <- mvrnorm(nTotal, mu = muCtrl, Sigma = B)
collaboRATE_intv <- mvrnorm(nTotal, mu = muIntv, Sigma = B)

Patient_Assignments <- subset(Setup6, Time2==0, select=c(Subject, Site_Assign))

collaboRATE_setup <- data.frame(matrix(ncol = nTimePts, nrow = nTotal))
for(b in 1:nTimePts){
for(a in 1:nTotal){
collaboRATE_setup[a,b] <- ifelse(Patient_Assignments[a,2]==0, collaboRATE_ctrl[a,b], collaboRATE_intv[a,b])
}
}

collaboRATE_wide=cbind(Patient_Assignments, collaboRATE_setup)

collaboRATE_long <- reshape(collaboRATE_wide, 
  varying = c("X2", "X3", "X4"), 
  v.names = "collaboRATE",
  timevar = "Time", 
  times = c("Time1", "Time2", "Time3"), 
  direction = "long")

collaboRATE_long_sort<-as.data.frame(collaboRATE_long[order(collaboRATE_long $Subject),])

Setup7 <- subset(Setup3, Subject>0, select=c(Subject, SurvivalTime, Site_Effect))

Setup8 = as.data.frame(merge(x=collaboRATE_long_sort, y=Setup7, by="Subject", all=TRUE))

Setup8$collaboRATE_bl <- Setup8$X1 + Setup8$Site_Effect
Setup8$collaboRATE_fu <- Setup8$collaboRATE + Setup8$Site_Effect
Setup8$collaboRATE_fu <- ifelse(Setup8$Time=="Time3" & Setup8$SurvivalTime<=180, NA, Setup8$collaboRATE_fu)
Setup8$collaboRATE_fu <- ifelse(Setup8$Time=="Time2" & Setup8$SurvivalTime<=90, NA, Setup8$collaboRATE_fu)
Setup8$collaboRATE_fu <- ifelse(Setup8$Time=="Time1" & Setup8$SurvivalTime<=30, NA, Setup8$collaboRATE_fu)


#####################################################################
#Final Dataset: Subject, Site, Site_Assign, collaboRATE_bl, Time, collaboRATE_fu

FinalSetup1 <- subset(Setup1, Subject>0, select=c(Subject, Site))

FinalSetup1$SiteAsFactor <- as.factor(FinalSetup1$Site)

FinalSetup2 <- subset(Setup3, Subject>0, select=c(Subject, Site_Assign))

FinalSetup3 = merge(x=FinalSetup1, y=FinalSetup2, by="Subject", all=TRUE)

FinalSetup4 <- subset(Setup8, Time=="Time1")
FinalSetup5 <- subset(FinalSetup4, Subject>0, select=c(Subject, collaboRATE_bl))

FinalSetup6 = merge(x=FinalSetup3 , y=FinalSetup5, by="Subject", all=TRUE)

FinalSetup7 <- subset(Setup8, Subject>0, select=c(Subject, Time, collaboRATE_fu))

FinalSetup8 = merge(x=FinalSetup6, y=FinalSetup7, by="Subject", all=TRUE)

sampletrial=na.omit(FinalSetup8)

  lm <- lm(collaboRATE_fu~Site,data=sampletrial)
  lm2=anova(lm)
  ICC[i]=(lm2[1,2])/(lm2[1,2]+lm2[2,2])

  model = lmer(collaboRATE_fu ~ Site_Assign + collaboRATE_bl + (1|Subject), data = sampletrial, REML=TRUE)
  b_1[i]=(summary(model)$coefficients[2])
  lcl_1[i]=(summary(model)$coefficients[2]-1.96*summary(model)$coefficients[5])  
  ucl_1[i]=(summary(model)$coefficients[2]+1.96*summary(model)$coefficients[5])
  pval_1[i]=summary(model)$coefficients[14]
  succ_1[i]=ifelse(lcl_1[i]>0, 1, 0)

  model2 = lmer(collaboRATE_fu ~ Site_Assign + SiteAsFactor + collaboRATE_bl + (1|Subject), data = sampletrial, REML=TRUE)
  b_2[i]=(summary(model2)$coefficients[2])
  lcl_2[i]=(summary(model2)$coefficients[2]-1.96*summary(model2)$coefficients[11])
  ucl_2[i]=(summary(model2)$coefficients[2]+1.96*summary(model2)$coefficients[11])
  pval_2[i]=summary(model2)$coefficients[38]
  succ_2[i]=ifelse(lcl_2[i]>0, 1, 0)

  model3 = lmer(collaboRATE_fu ~ Site_Assign  + collaboRATE_bl + (1|Site/Subject), data = sampletrial, REML=TRUE)
  b_3[i]=(summary(model3)$coefficients[2])
  lcl_3[i]=(summary(model3)$coefficients[2]-1.96*summary(model3)$coefficients[5])
  ucl_3[i]=(summary(model3)$coefficients[2]+1.96*summary(model3)$coefficients[5])
  pval_3[i]=summary(model3)$coefficients[14]
  succ_3[i]=ifelse(lcl_3[i]>0, 1, 0)
}

#return(sampletrial)

#trialresults=as.data.frame(cbind(trial, sd_site, ICC, b_1, lcl_1, ucl_1, pval_1, succ_1))
#return(trialresults)

#trialresults=as.data.frame(cbind(trial, sd_site, ICC, b_1, lcl_1, ucl_1, pval_1, succ_1, b_2, lcl_2, ucl_2, pval_2, succ_2))
#return(trialresults)

trialresults=as.data.frame(cbind(trial, sd_site, ICC, b_1, lcl_1, ucl_1, pval_1, succ_1, b_2, lcl_2, ucl_2, pval_2, succ_2, b_3, lcl_3, ucl_3, pval_3, succ_3))
return(trialresults)

}

#sampletrial=TrialSimulation(seed=1, nSims=1, nSites=8, nPerSite=75, nTimePts=4, site_var=3,
#muCtrl=c(100, 100, 100, 100), 
#muIntv=c(100, 103, 103, 103),
#stdevs=c(6, 6, 6, 6), 
#mvncorr=0.75, weibull_param=1700)
#head(sampletrial, n=10)

#model2 = lmer(collaboRATE_fu ~ Site_Assign + SiteAsFactor + collaboRATE_bl + (1|Subject), data = sampletrial, REML=TRUE)
#summary(model2)
#b_2=(summary(model2)$coefficients[2])
#b_2
#lcl_2=(summary(model2)$coefficients[2]-1.96*summary(model2)$coefficients[11])
#lcl_2
#ucl_2=(summary(model2)$coefficients[2]+1.96*summary(model2)$coefficients[11])
#ucl_2  
#summary(model2)
#pval_2=summary(model2)$coefficients[38]
#pval_2

trialresults1=TrialSimulation(seed=1, nSims=10, nSites=6, nPerSite=75, nTimePts=4, site_var=6,
muCtrl=c(100, 100, 100, 100), 
muIntv=c(100, 103, 103, 103),
stdevs=c(6, 6, 6, 6), 
mvncorr=0.75, weibull_param=1700)
head(trialresults1, n=10)
summary(trialresults1$sd_site)
summary(trialresults1$ICC)
table(trialresults1$succ_1)
table(trialresults1$succ_2)


