ctrl_values<-as.numeric(c(1,4,0,2,0,7,9,8,4,12,0,8,6,0,0,7,1,0,1,0,9,1,14,11,0,0,21,8,0,0,9,0,16,14,5,17,9,23,3,1,7,3,2,9,0,7,0,2,0,18,3,3,24,6,19,2,1))

trial_simulation <- function(nSims, nPerGroup, control_values, trt_effect, seed) {

set.seed(seed)

trial <- numeric(nSims)
p_ttest <- numeric(nSims)
success_t <- numeric(nSims)
p_wilcoxon <- numeric(nSims)
success_w <- numeric(nSims)

for(i in 1:nSims){

trial[i]=i

ctrl_dist<-control_values
trt_dist<-ctrl_dist*(1-trt_effect)

#Create Sample Trial Data
control <- sample(ctrl_dist,nPerGroup,replace=TRUE)
treat <- sample(trt_dist,nPerGroup,replace=TRUE)

#Analyses
t <- t.test(control, treat)
p_ttest[i] <- round(t$p.value,4)
success_t[i] <- ifelse(mean(treat) < mean(control) & p_ttest[i] < 0.05, 1, 0)
w <- wilcox.test(control, treat)
p_wilcoxon[i] <- round(w$p.value,4)
success_w[i] <- ifelse(mean(treat) < mean(control) & p_wilcoxon[i] < 0.05, 1, 0)
}

trialresults=as.data.frame(cbind(trial, p_ttest, success_t, p_wilcoxon, success_w))

return(trialresults)

}

trialresults1=trial_simulation(nSims=1000, nPerGroup=50, control_values=ctrl_values, trt_effect=0.5, seed=1)
trialresults2=trial_simulation(nSims=1000, nPerGroup=75, control_values=ctrl_values, trt_effect=0.5, seed=1)
trialresults3=trial_simulation(nSims=1000, nPerGroup=100, control_values=ctrl_values, trt_effect=0.5, seed=1)
trialresults4=trial_simulation(nSims=1000, nPerGroup=125, control_values=ctrl_values, trt_effect=0.5, seed=1)
trialresults5=trial_simulation(nSims=1000, nPerGroup=150, control_values=ctrl_values, trt_effect=0.5, seed=1)

table(trialresults1$success_t)
table(trialresults2$success_t)
table(trialresults3$success_t)
table(trialresults4$success_t)
table(trialresults5$success_t)

table(trialresults1$success_w)
table(trialresults2$success_w)
table(trialresults3$success_w)
table(trialresults4$success_w)
table(trialresults5$success_w)
