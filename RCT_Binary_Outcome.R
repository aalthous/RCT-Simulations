# This code will mimic a 2-group parallel-arm randomized trial using 1:1 allocation of patients to treatment 1 versus treatment 2
# For this example, we will use a binary outcome of "death"
# Patients receiving treatment 1 will have 40% probability of death
# Patients receiving treatment 2 will have 30% probability of death
# Analysis will be performed using a logistic regression model with no covariate adjustment
# We will run 1000 simulated RCT's and report the odds ratio, 95% confidence interval, and p-value for each simulated trial
# The "true" treatment effect for a treatment that reduces probability of outcome from 40% to 30% is about OR = 0.642
# The power of a trial with N=1000 patients and exactly 1:1 allocation under these assumptions is about 91-92%

# Trial Design Parameters
nPatients <- 1000 # here is where you specify the number of patients you want included in each RCT 
death1 <- 0.4 # here is where you specify the event rate for patients receiving 'treatment 1' in these trials
death2 <- 0.3 # here is where you specify the event rate for patients receiving 'treatment 2' in these trials

# Simulation Parameters
nSims <- 1000 # here is where you specify the number of trials that you want to simulate
trialnum <- numeric(nSims) # this creates an empty vector that we'll populate counting upwards from 1 to nSims
or <-numeric(nSims) # this creates an empty vector that we'll populate with the odds ratio estimate for each simulated trial
lcl <-numeric(nSims) # this creates an empty vector that we'll populate with the lower limit of 95% CI for the odds ratio for each simulated trial
ucl <-numeric(nSims) # this creates an empty vector that we'll populate with the upper limit of 95% CI for the odds ratio for each simulated trial
pvalue <-numeric(nSims) # this creates an empty vector that we'll populate with the p-value for each simulated trial
success <-numeric(nSims) # this creates an empty vector that we'll populate with an indicator (1 for "success" / 0 for "failure") for each simulated trial

set.seed(1) # this sets the random seed for your results to be reproducible

for(i in 1:nSims){

pid=seq(1, by=1, len=nPatients) # this creates a sequential list of "pid" from 1 to nPatients which may be useful if you want to perform 'interim analysis' later
treatment=rep(1:2, nPatients/2) # this creates a vector of "treatment allocations" which is actually just a sequence alternating between 1 and 2

# worth noting: this allocation sequence should not be used in a real RCT, but for the purpose of these simulations it will work fine.  
# There are no real patients or clinicians created in these simulations, and therefore no worry about someone guessing the next treatment assignment.
# If you prefer that your simulations actually use “randomized” allocation, you can do this instead:
# treatment=1+rbinom(nPatients, 1, 0.5) # this randomly assigns each new patient to receive treatment 1 or 2 with 50% probability each time
# The reason I prefer the first of the two for simulation is that it maintains even allocation in the number of patients receiving each treatment 
# (of course, with a wee bit more work one can actually create blocked randomization sequence, but I’m trying to keep this simple for newbies)
# (for those interested in going one step further, the "blockrand" package can be used to generate this, may include in future posts)
# The "simple randomization" example will have *slightly* lower power due to the allowance for an imbalanced number of patients; 
# Using "exactly-equal-allocation" means we will be *slightly* over-estimating the trial power by assuming exactly equal allocation
# when stratified and/or blocked randomization could allow slightly unequal allocations to occur, e.g. "498 vs 502" patients
# Constraining to "exactly equal" is close enough in practice to results with blocked randomization that it's my preference
# Also worth noting, most people doing conventional power calculations (without simulation) assume exactly equal allocation
# without accounting for the slight imbalances that may occur in truly 'random' allocation sequences

deathprob <- numeric(nPatients) # this creates an empty vector which we will use to assign death probability for each patient
deathprob[treatment==1]=death1 # this assigns the probability of death for patients receiving 'treatment 1' to be 'death1'
deathprob[treatment==2]=death2 # this assigns the probability of death for patients receiving 'treatment 2' to be 'death2'
death=rbinom(nPatients, 1, deathprob) # this simulates each patient's outcome as a random draw from binomial distribution with probabilities assigned above
trialdata=data.frame(cbind(pid, treatment, death)) # this creates a data frame with pid, treatment allocation, and death outcome

model <- glm(death ~ treatment, family=binomial(link='logit'), data=trialdata) # this runs a logistic regression model on each trial's simulated data
trialnum[i]=i # this simply tells us which simulation each row of our results came from (counting upward from 1 to nSims)
or[i]=round(exp(summary(model)$coefficients[2]), digits=2) # this saves the odds ratio for the treatment effect from the logistic regression model
lcl[i]=round(exp(summary(model)$coefficients[2] - 1.96 * summary(model)$coefficients[4]), digits=2) # this saves the lower limit of a 95% CI for the treatment effect OR
ucl[i]=round(exp(summary(model)$coefficients[2] + 1.96 * summary(model)$coefficients[4]), digits=2) # this saves the upper limit of a 95% CI for the treatment effect OR
pvalue[i]=round(summary(model)$coefficients[8], digits=4) # this saves the p-value for the treatment effect
success[i]=ifelse(or[i]<1 & pvalue[i]<0.05, 1, 0) # this creates a flag for whether the trial was a "success" - if effect favored "treatment 2" and p<0.05 for treatment effect

}

simulation_results <- data.frame(cbind(trialnum, or, lcl, ucl, pvalue, success)) # this puts all of the saved results into one data frame for easy viewing
head(simulation_results , n=10) # this lets you take a look at the first 10 simulation results to confirm all seems to make sense
table(success) # this provides a table of the number of simulated trials which concluded that there was a treatment effect


