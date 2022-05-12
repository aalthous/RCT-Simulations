# RCT-Simulation-v1
# This code will mimic a 2-group parallel-arm randomized trial using 1:1 allocation of patients to treatment 1 versus treatment 2
# For this example, we will use a binary outcome of "death"
# Patients receiving treatment 1 will have 40% probability of death
# Patients receiving treatment 2 will have 30% probability of death
# Analysis will be performed using a logistic regression model for interim and final analysis
# We will run 1000 simulated RCT's and report the odds ratio, 95% confidence interval, and p-value for each simulated trial at interim and final analysis
# The "true" treatment effect for a treatment that reduces probability of outcome from 40% to 30% is about OR = 0.642
# The power of a trial with N=1000 patients and exactly 1:1 allocation under these assumptions is about 91-92%
# In this example, we are now adding code that mimics performing one interim analysis with an efficacy stopping rule
# The interim strategy in the code shown here mimics an O'Brien-Fleming rule with one interim analysis that takes place at 50% observed information
# The power of a trial with N=1000 patients and exactly 1:1 allocation under these assumptions with this design is about 91-92%

# Trial Design Parameters
nPatients <- 1000 # here is where you specify the planned max number of patients you want included in each RCT 
nInterim <- 500 # here is where you specify the number of patients where you want the interim anaysis to occur
# NOTE: the alpha thresholds specified here are set for the O'Brien Fleming rule with an interim occurring at 50% observed information
#   if you want to use a different time for the interim analysis (e.g. 75% rather than 50%) need to adjust the alpha thresholds below
interim_efficacy_threshold <- 0.0054 # here is where you specify the p-value for stopping at interim analysis 
final_efficacy_threshold <- 0.0492 # here is where you specify the p-value for declaring success at final if trial continues past the interim analysis
# NOTE: 0.0054 and 0.0492 are the thresholds for an O'Brien Fleming with one interim analysis that takes place at 50% of the observed outcome data
# If you want to fiddle with different stopping thresholds, I do suggest reading up on different approaches to alpha spending, but the simulations here
#    will be interesting for you to play with if you want to see how operating characteristics (type 1 & type 2 error) change with different thresholds
#    I will explore this a bit more in a subsequent post that uses the "rpact" package to compute the efficacy stopping thresholds
death1 <- 0.4 # here is where you specify the event rate for patients receiving 'treatment 1' in these trials
death2 <- 0.3 # here is where you specify the event rate for patients receiving 'treatment 2' in these trials
# NOTE: if you want to estimate "type 1 error" under different stopping rules, make death = in the two treatment arms (e.g. no treatment effect)
#    this can be useful to get across that type 1 error increases if you take multiple looks at data using p<0.05 without appropriately accounting 
#    for the multiple looks at the data

# Simulation Parameters
nSims <- 1000 # here is where you specify the number of trials that you want to simulate
trialnum <- numeric(nSims) # this creates an empty vector that we'll populate counting upwards from 1 to nSims
or_interim <-numeric(nSims) # this creates an empty vector that we'll populate with the odds ratio estimate for each simulated trial
lcl_interim <-numeric(nSims) # this creates an empty vector that we'll populate with the lower limit of 95% CI for the odds ratio for each simulated trial
ucl_interim <-numeric(nSims) # this creates an empty vector that we'll populate with the upper limit of 95% CI for the odds ratio for each simulated trial
pvalue_interim <-numeric(nSims) # this creates an empty vector that we'll populate with the p-value for each simulated trial
success_interim <-numeric(nSims) # this creates an empty vector that we'll populate with an indicator (1 for "success" / 0 for "failure") for each simulated trial
or_final <-numeric(nSims) # this creates an empty vector that we'll populate with the odds ratio estimate for each simulated trial
lcl_final <-numeric(nSims) # this creates an empty vector that we'll populate with the lower limit of 95% CI for the odds ratio for each simulated trial
ucl_final <-numeric(nSims) # this creates an empty vector that we'll populate with the upper limit of 95% CI for the odds ratio for each simulated trial
pvalue_final <-numeric(nSims) # this creates an empty vector that we'll populate with the p-value for each simulated trial
success_final <-numeric(nSims) # this creates an empty vector that we'll populate with an indicator (1 for "success" / 0 for "failure") for each simulated trial
overall_success <- numeric(nSims) # this creates an empty vector that we'll populate with an indicator (1 for "success" / 0 for "failure") for each simulated trial

set.seed(1) # this sets the random seed for your results to be reproducible

for(i in 1:nSims){

pid=seq(1, by=1, len=nPatients) # this creates a sequential list of "pid" from 1 to nPatients which may be useful if you want to perform 'interim analysis' later
treatment=rep(1:2, nPatients/2) # this creates a vector of "treatment allocations" which is actually just a sequence alternating between 1 and 2

# worth noting: this allocation sequence should not be used in a real RCT, but for the purpose of these simulations it will work fine.  
# There are no real patients or clinicians created in these simulations, and therefore no worry about someone guessing the next treatment assignment.
# If you prefer that your simulations actually use “randomized” allocation, you can do this instead:
# treatment=1+rbinom(nPatients, 1, 0.5) # this randomly assigns each new patient to receive treatment 1 or 2 with 50% probability each time
# The reason I prefer the first of the two for simulation is that it maintains even allocation in the number of patients receiving each treatment 
# (of course, with a wee bit more work one can actually create blocked randomization sequence, but I’m trying to keep this thread simple for newbies)
# (for those interested in going one step further, the "blockrand" package can be used to generate this, may include in future post)
# The "simple randomization" example will have slightly lower power due to the allowance for an imbalanced number of patients; 
# Using "exactly-equal-allocation" means we will be slightly over-estimating the trial power by assuming exactly equal allocation
# when stratified and/or blocked randomization could allow slightly unequal allocations to occur, e.g. "498 vs 502" patients
# Constraining to "exactly equal" is close enough in practice to results with blocked randomization that it's my preference

deathprob <- numeric(nPatients) # this creates an empty vector which we will use to assign death probability for each patient
deathprob[treatment==1]=death1 # this assigns the probability of death for patients receiving 'treatment 1' to be 'death1'
deathprob[treatment==2]=death2 # this assigns the probability of death for patients receiving 'treatment 2' to be 'death2'
death=rbinom(nPatients, 1, deathprob) # this simulates each patient's outcome as a random draw from binomial distribution with probabilities assigned above
trialdata=data.frame(cbind(pid, treatment, death)) # this creates a data frame with pid, treatment allocation, and death outcome

interimdata<- subset(trialdata, pid<=nInterim) # this selects just the patients with outcomes known at time of interim analysis

trialnum[i]=i # this simply tells us which simulation each row of our results came from (counting upward from 1 to nSims)

model1 <- glm(death ~ treatment, family=binomial(link='logit'), data=interimdata) # this runs a logistic regression model on the interim analysis set
or_interim[i]=round(exp(summary(model1)$coefficients[2]), digits=2) # this saves the odds ratio for the treatment effect from the logistic regression model
lcl_interim[i]=round(exp(summary(model1)$coefficients[2] - 1.96 * summary(model1)$coefficients[4]), digits=2) # this saves the lower limit of a 95% CI for the treatment effect OR
ucl_interim[i]=round(exp(summary(model1)$coefficients[2] + 1.96 * summary(model1)$coefficients[4]), digits=2) # this saves the upper limit of a 95% CI for the treatment effect OR
pvalue_interim[i]=round(summary(model1)$coefficients[8], digits=4) # this saves the p-value for the treatment effect at interim analysis
success_interim[i]=ifelse(or_interim[i]<1 & pvalue_interim[i]<interim_efficacy_threshold, 1, 0) # this creates a flag for whether the trial was a "success" - if effect favored "treatment 2" and p<0.05 for treatment effect

model2 <- glm(death ~ treatment, family=binomial(link='logit'), data=trialdata) # this runs a logistic regression model on the full analysis set
or_final[i]=round(exp(summary(model2)$coefficients[2]), digits=2) # this saves the odds ratio for the treatment effect from the logistic regression model
lcl_final[i]=round(exp(summary(model2)$coefficients[2] - 1.96 * summary(model2)$coefficients[4]), digits=2) # this saves the lower limit of a 95% CI for the treatment effect OR
ucl_final[i]=round(exp(summary(model2)$coefficients[2] + 1.96 * summary(model2)$coefficients[4]), digits=2) # this saves the upper limit of a 95% CI for the treatment effect OR
pvalue_final[i]=round(summary(model2)$coefficients[8], digits=4) # this saves the p-value for the treatment effect at final analysis
success_final[i]=ifelse(or_final[i]<1 & pvalue_final[i]<final_efficacy_threshold, 1, 0) # this creates a flag for whether the trial was a "success" - if effect favored "treatment 2" and p<0.05 for treatment effect

overall_success[i]=ifelse(success_interim[i] == 1, 1, success_final[i])

}

simulation_results <- data.frame(cbind(trialnum, or_interim, lcl_interim, ucl_interim, pvalue_interim, success_interim,
or_final, lcl_final, ucl_final, pvalue_final, success_final, overall_success)) # this puts all of the saved results into one data frame for easy viewing
head(simulation_results , n=10) # this lets you take a look at the first 10 simulation results to confirm all seems to make sense
table(success_interim) # this provides a table of the number of simulated trials which concluded success at the interim analysis
table(overall_success) # this provides a table of the number of simulated trials which concluded that there was a treatment effect
