# Trial Design Parameters - Part 1
# Here we will specify the basics: total N patients to enroll, and death rate for each treatment arm
nPatients <- 1000 # here is where you specify the planned max number of patients you want included in each RCT 
death1 <- 0.4 # here is where you specify the event rate for patients receiving 'treatment 1' in these trials
death2 <- 0.3 # here is where you specify the event rate for patients receiving 'treatment 2' in these trials
# NOTE: if you want to confirm "type 1 error" under different stopping rules, make death = in the two treatment arms (e.g. no treatment effect)
# I have set this one up to test the power for a treatment that would reduce mortality from 40% in control group (1) to 30% in treatment group (2)

# Trial Design Parameters - Part 2
# Here we will define the interim analysis strategy and stopping rules
# For this trial we will include provisions for efficacy stopping only (no pre-specified futility stopping)
# We will use the rpact package to compute the stopping/success thresholds at the interim and final analysis 
# install.packages("rpact")
# library(rpact)
nLooks<-4 # here is where you put the number of looks that will take place (INCLUDING the final analysis)
analyses_scheduled<-(c(0.25, 0.50, 0.75, 1))
efficacy_thresholds<-numeric(nLooks)

design <- getDesignGroupSequential(sided=1, alpha=0.05, informationRates=analyses_scheduled, typeOfDesign = "asOF")
for(j in 1:nLooks){
efficacy_thresholds[j] = design$stageLevels[j]
}
analyses_scheduled
analyses_nPatients <- analyses_scheduled*nPatients
analyses_nPatients
efficacy_thresholds

# Simulation Parameters
nSims <- 1000
trialnum <- numeric(nSims)

or <- data.frame(matrix(ncol = nLooks, nrow = nSims))
lcl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
ucl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
pvalue <- data.frame(matrix(ncol = nLooks, nrow = nSims))
success <- data.frame(matrix(ncol = nLooks, nrow = nSims))

#provide column names
colnames(or) <- c('or_1', 'or_2', 'or_3', 'or_4')
colnames(lcl) <- c('lcl_1', 'lcl_2', 'lcl_3', 'lcl_4')
colnames(ucl) <- c('ucl_1', 'ucl_2', 'ucl_3', 'ucl_4')
colnames(pvalue) <- c('pvalue_1', 'pvalue_2', 'pvalue_3', 'pvalue_4')
colnames(success) <- c('success_1', 'success_2', 'success_3', 'success_4')

overall_success <- numeric(nSims)

set.seed(1) # this sets the random seed for your results to be reproducible

for(i in 1:nSims){

trialnum[i]=i

pid=seq(1, by=1, len=nPatients)
treatment=rep(1:2, nPatients/2)
deathprob <- numeric(nPatients)
deathprob[treatment==1]=death1
deathprob[treatment==2]=death2
death=rbinom(nPatients, 1, deathprob)
trialdata=data.frame(cbind(pid, treatment, death))

for(j in 1:nLooks){
analysisdata <- subset(trialdata, pid<=analyses_nPatients[j])
model <- glm(death ~ treatment, family=binomial(link='logit'), data=analysisdata)
or[i,j]=exp(summary(model)$coefficients[2])
lcl[i,j]=exp(summary(model)$coefficients[2] - 1.96 * summary(model)$coefficients[4])
ucl[i,j]=exp(summary(model)$coefficients[2] + 1.96 * summary(model)$coefficients[4]) 
pvalue[i,j]=summary(model)$coefficients[8]
success[i,j]=ifelse(or[i,j]<1 & pvalue[i,j]<efficacy_thresholds[j], 1, 0)
}

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
or[1], lcl[1], ucl[1], pvalue[1], success[1],
or[2], lcl[2], ucl[2], pvalue[2], success[2],
or[3], lcl[3], ucl[3], pvalue[3], success[3],
or[4], lcl[4], ucl[4], pvalue[4], success[4],
overall_success))
head(simulation_results, n=10)

table(overall_success)
table(simulation_results$success_1, overall_success)
table(simulation_results$success_2, overall_success)
table(simulation_results$success_3, overall_success)
table(simulation_results$success_4, overall_success)


