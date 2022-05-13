# This code will mimic a 2-group parallel-arm randomized trial using 1:1 allocation of patients to treatment 0 ("control") versus treatment 1 ("treatment")
# For this example, we will use a time to event outcome
# Patients receiving treatment 0 will have a survival time drawn from Weibull distribution with shape parameter 0.5 and scale parameter 1200
# This corresponds to approximately 70.6% of patients experiencing the event by t=1800 (5 years) from that patient's date of randomization
# Patients receiving treatment 1 will have a survival time drawn from Weibull distribution with shape parameter 0.5 and scale parameter 2000
# This corresponds to approximately 61.3% of patients experiencing the event by t=1800 (5 years) from that patient's date of randomization
# There are also steps baked in below that simulate pace of recruitment and censoring of patients "still alive" at designated end of trial 
# In this example, I am *ending* the trial at t=1800 (5 years) from the start of the trial, which means the longest possible follow-up time
#      is 5 years, and therefore patients enrolled later in the trial will have less available follow-up time (mimicking 'real life' trial)
# Analysis will be performed using a Cox proportional-hazards regression model 
# We will run 1000 simulated RCT's and report the hazard ratio, 95% confidence interval, and p-value for each simulated trial
# The "true" treatment effect for a treatment with these survival rates is about HR = 0.xx
# The power of a trial with N=1000 patients and exactly 1:1 allocation under these assumptions is about xx%
# IMPORTANT: if you are going to tinker with this, you need to be mindful to adjust all of the "times" appropriately to mimic your trial
# IMPORTANT: related to above, my times are scaled in "days" in this example; if you are thinking "months" or "years" need to be careful
# IMPORTANT: finally, understanding what the shape and scale parameters do for simulation of survival times is a bit beyond intro scope
# IMPORTANT: for beginners, I think the best thing to do is use sample code below to view Kaplan-Meier curves under different values of
#            the shape and scale parameters to understand the implications of different choices in these values

# Before we get into the full simulation code, here it can be very helpful to create one toy dataset to visualize
# what the survival curves look like from your desired sample, so before starting the full simulation code (below) 
# we'll make one sample dataset and plot some Kaplan-Meier curves

set.seed(1)
library(survival) # Sorry, this time you have to load at least one package; no getting around it
nPerGroup <- 500 # Number of patients per treatment group
weibull_shape0 <- 0.5 # Shape parameter for Weibull in treatment group 0
weibull_shape1 <- 0.5 # Shape parameter for Weibull in treatment group 1
weibull_scale0 <- 1200 # Scale parameter for Weibull in treatment group 0
weibull_scale1 <- 2000  # Scale parameter for Weibull in treatment group 1
final_analysis_date <- 1800 # Final "analysis date" (for this example, t=1800 means dataset is frozen 5 years after trial begins)

#Create One Sample Trial Dataset
pid=seq(1, by=1, len=nPerGroup*2) # this creates a sequential list of "pid" from 1 to nPatients which may be useful if you want to perform 'interim analysis' later
treatment=rep(0:1, nPerGroup) # this creates a vector of "treatment allocations" which is actually just a sequence alternating between 0 and 1
recruitdate=floor((pid+1)/2) # create recruitment date; this mimics 2 patients recruited daily (a trial of n=500 per group will reach target in 500 days)
survivaltime <- numeric(nPerGroup*2) # create empty vector which will hold survival times
survivaltime[treatment==0] <- rweibull(nPerGroup, weibull_shape0, weibull_scale0) # simulated survival times from Weibull distribution for treatment 0
survivaltime[treatment==1] <- rweibull(nPerGroup, weibull_shape1, weibull_scale1) # simulated survival times from Weibull distribution for treatment 1
sampletrial=data.frame(cbind(pid, treatment, recruitdate, survivaltime)) # merge patient id, treatment allocation, recruit date, and simulated survival time
head(sampletrial, n=10) # take a look at this to see what we've created so far

# Create Indicator Variables For Final Analysis
# All patients will have their status and survival time computed as of "final_analysis_date" specified above

sampletrial$deathdate=sampletrial$recruitdate+sampletrial$survivaltime # this creates a "death date" which is just the survival time added to recruitment date

sampletrial$death[sampletrial$deathdate<final_analysis_date]=1 # if the death date is less than final analysis date, death=1
sampletrial$death[sampletrial$deathdate>=final_analysis_date]=0 # if the death date is greater than final analysis date, death=1
head(sampletrial, n=10) # take a look at this to see what we've created so far

sampletrial$deathtime<- 
ifelse(sampletrial$death == 1,
sampletrial$deathdate-sampletrial$recruitdate,
final_analysis_date-sampletrial$recruitdate) #this creates a "deathtime" variable which represents survival time as of "final analysis date"

head(sampletrial, n=10) # take a look at this to see what we've created so far

#Kaplan-Meier Curve for Sample Trial Data

KM_final <- survfit(Surv(deathtime, death) ~ treatment, type="kaplan-meier", conf.type=c("none"), data=sampletrial)
plot(KM_final, main=expression(paste("Kaplan-Meier Estimate")), xlab="Time (Years)", ylab="Survival", lwd=2, col=1:2)

coxmodel=coxph(Surv(deathtime, death) ~ treatment, data=sampletrial)
summary(coxmodel)

# Okay, do those curves look reasonably like what you expected?  If so, then proceed with code below
# If not, tinker with the shape and scale parameters of the Weibull until the sample K-M curves look
# like the curves you imagine represent the outcomes you want to use for your simulations

# Now we'll take that code, wrap it into a for loop, and run a bunch of those simulated trials
# We'll skip the Kaplan-Meier curve in the for loop, since that was just for illustration

# Trial Design Parameters
nPerGroup <- 500
weibull_shape0 <- 0.5
weibull_shape1 <- 0.5
weibull_scale0 <- 1200
weibull_scale1 <- 2000
final_analysis_date <- 1800 # Final "analysis date" (for this example, t=1800 means dataset is frozen 5 years after trial begins)

# Simulation Parameters
nSims <- 1000 # here is where you specify the number of trials that you want to simulate
trialnum <- numeric(nSims) # this creates an empty vector that we'll populate counting upwards from 1 to nSims
hr <-numeric(nSims) # this creates an empty vector that we'll populate with the hazard ratio estimate for each simulated trial
lcl <-numeric(nSims) # this creates an empty vector that we'll populate with the lower limit of 95% CI for the hazard ratio for each simulated trial
ucl <-numeric(nSims) # this creates an empty vector that we'll populate with the upper limit of 95% CI for the hazard ratio for each simulated trial
pvalue <-numeric(nSims) # this creates an empty vector that we'll populate with the p-value for each simulated trial
success <-numeric(nSims) # this creates an empty vector that we'll populate with an indicator (1 for "success" / 0 for "failure") for each simulated trial

set.seed(1) # this sets the random seed for your results to be reproducible

for(i in 1:nSims){

pid=seq(1, by=1, len=nPerGroup*2) # this creates a sequential list of "pid" from 1 to nPatients which may be useful if you want to perform 'interim analysis' later
treatment=rep(0:1, nPerGroup) # this creates a vector of "treatment allocations" which is actually just a sequence alternating between 0 and 1

# worth noting: this allocation sequence should not be used in a real RCT, but for the purpose of these simulations it will work fine.  
# There are no real patients or clinicians created in these simulations, and therefore no worry about someone guessing the next treatment assignment.
# If you prefer that your simulations actually use “randomized” allocation, you can do this instead:
# treatment=rbinom(nPatients, 1, 0.5) # this randomly assigns each new patient to receive treatment 0 or 1 with 50% probability each time
# The reason I prefer the first of the two for simulation is that it maintains even allocation in the number of patients receiving each treatment 
# (of course, with a wee bit more work one can actually create blocked randomization sequence, but I’m trying to keep this simple for newbies)
# (for those interested in going one step further, the "blockrand" package can be used to generate this, may include in future posts)
# The "simple randomization" example will have *slightly* lower power due to the allowance for an imbalanced number of patients; 
# Using "exactly-equal-allocation" means we will be *slightly* over-estimating the trial power by assuming exactly equal allocation
# when stratified and/or blocked randomization could allow slightly unequal allocations to occur, e.g. "498 vs 502" patients
# Constraining to "exactly equal" is close enough in practice to results with blocked randomization that it's my preference
# Also worth noting, most people doing conventional power calculations (without simulation) assume exactly equal allocation
# without accounting for the slight imbalances that may occur in truly 'random' allocation sequences

recruitdate=floor((pid+1)/2) # create recruitment date; this mimics 2 patients recruited daily (a trial of n=500 per group will reach target in 500 days)
survivaltime <- numeric(nPerGroup*2) # create empty vector which will hold survival times
survivaltime[treatment==0] <- rweibull(nPerGroup, weibull_shape0, weibull_scale0) # simulated survival times from Weibull distribution for treatment 0
survivaltime[treatment==1] <- rweibull(nPerGroup, weibull_shape1, weibull_scale1) # simulated survival times from Weibull distribution for treatment 1
trialdata=data.frame(cbind(pid, treatment, recruitdate, survivaltime)) # merge patient id, treatment allocation, recruit date, and simulated survival time

# Create Indicator Variables For Final Analysis
# All patients will have their status and survival time computed as of "final_analysis_date" specified above

trialdata$deathdate=trialdata$recruitdate+trialdata$survivaltime # this creates a "death date" which is just the survival time added to recruitment date

trialdata$death[trialdata$deathdate<final_analysis_date]=1 # if the death date is less than final analysis date, death=1
trialdata$death[trialdata$deathdate>=final_analysis_date]=0 # if the death date is greater than final analysis date, death=1

trialdata$deathtime<- 
ifelse(trialdata$death == 1,
trialdata$deathdate-trialdata$recruitdate,
final_analysis_date-trialdata$recruitdate) #this creates a "deathtime" variable which represents survival time as of "final analysis date"

coxmodel=coxph(Surv(deathtime, death) ~ treatment, data=trialdata)
trialnum[i]=i # this simply tells us which simulation each row of our results came from (counting upward from 1 to nSims)
hr[i]=round(exp(summary(coxmodel)$coefficients[1]), digits=2) # this saves the hazard ratio for the treatment effect from the Cox proportional-hazards model
lcl[i]=round(exp(summary(coxmodel)$coefficients[1]-1.96*summary(coxmodel)$coefficients[3]), digits=2)  # this saves the lower limit of a 95% CI for the treatment effect HR
ucl[i]=round(exp(summary(coxmodel)$coefficients[1]+1.96*summary(coxmodel)$coefficients[3]), digits=2) # this saves the upper limit of a 95% CI for the treatment effect HR
pvalue[i]=round(summary(coxmodel)$coefficients[5], digits=4) # this saves the p-value for the treatment effect
success[i]=ifelse(hr[i]<1 & pvalue[i]<0.05, 1, 0) # this creates a flag for whether the trial was a "success" - if effect favored "treatment 2" and p<0.05 for treatment effect

}

simulation_results <- data.frame(cbind(trialnum, hr, lcl, ucl, pvalue, success)) # this puts all of the saved results into one data frame for easy viewing
head(simulation_results , n=10) # this lets you take a look at the first 10 simulation results to confirm all seems to make sense
table(success) # this provides a table of the number of simulated trials which concluded that there was a treatment effect

