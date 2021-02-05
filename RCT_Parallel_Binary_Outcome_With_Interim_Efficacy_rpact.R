#' Simulation of a 2-group parallel-arm randomized trial with interim analyses using the rpact package
#' @author Andrew Althouse
#' @param nPatients here is where you specify the planned max number of patients you want included in each RCT
#' @param death1 here is where you specify the event rate for patients receiving 'treatment 1' in these trial
#' @param death2 # here is where you specify the event rate for patients receiving 'treatment 2' in these trials
#' @param nLooks # here is where you put the number of looks that will take place (INCLUDING the final analysis)
#' @param analyses_scheduled # schedule of interim analyses
#' @param sided # Whether the test is 1-sided or 2-sided
#' @param alpha # Specified alpha level, the default is 0.05
#' @param informationRates
#' @param typeOfDesign # The type of design. Type of design is one of the following: O'Brien & Fleming ("OF"), Pocock ("P"), Wang & Tsiatis Delta class ("WT"), Haybittle & Peto ("HP"), Optimum design within Wang & Tsiatis class ("WToptimum"), O'Brien & Fleming type alpha spending ("asOF"), Pocock type alpha spending ("asP"), Kim & DeMets alpha spending ("asKD"), Hwang, Shi & DeCani alpha spending ("asHSD"), user defined alpha spending ("asUser"), default is "OF".
#' @param nSims # Number of planned simulations
#' @param seed # Argument to set the seed for the simulations
#' @return
#' @export
#'
#' @examples
#'
rct_interim_sim <- function(nPatients = 1000, death1 = 0.4, death2 = 0.3, nLooks = 4, analyses_scheduled = c(0.25, 0.50, 0.75, 1),
                            sided = 1, alpha = 0.05, informationRates = analyses_scheduled, typeOfDesign = "asOF", nSims = 1000, seed = 1031) {

  # Trial Design Parameters - Part 1
  # Here we will specify the basics: total N patients to enroll, and death rate for each treatment arm
  seed <- seed
  nPatients <- nPatients # here is where you specify the planned max number of patients you want included in each RCT
  death1 <- death1 # here is where you specify the event rate for patients receiving 'treatment 1' in these trials
  death2 <- death2 # here is where you specify the event rate for patients receiving 'treatment 2' in these trials
  # NOTE: if you want to confirm "type 1 error" under different stopping rules, make death = in the two treatment arms (e.g. no treatment effect)
  # I have set this one up to test the power for a treatment that would reduce mortality from 40% in control group (1) to 30% in treatment group (2)

  # Trial Design Parameters - Part 2
  # Here we will define the interim analysis strategy and stopping rules
  # For this trial we will include provisions for efficacy stopping only (no pre-specified futility stopping)
  # We will use the rpact package to compute the stopping/success thresholds at the interim and final analysis
  # install.packages("rpact")
  library("rpact")
  nLooks <- nLooks # here is where you put the number of looks that will take place (INCLUDING the final analysis)

  sprintf("x%d", (1:nLooks))
  analyses_scheduled <- analyses_scheduled
  efficacy_thresholds <- numeric(nLooks)

  design <- getDesignGroupSequential(
    sided = sided, alpha = alpha,
    informationRates = analyses_scheduled, typeOfDesign = typeOfDesign)
  d2 <- getDesignGroupSequential(typeOfDesign = "P")
  d3 <- getDesignGroupSequential(typeOfDesign = "asP")
  d4 <- getDesignGroupSequential(typeOfDesign = "OF")

  designSet <- getDesignSet(
    designs = c(design, d2, d3, d4),
    variedParameters = "typeOfDesign")

  set.seed(seed)
  for (j in 1:nLooks) {
    efficacy_thresholds[j] <- design$stageLevels[j]
  }
  analyses_scheduled
  analyses_nPatients <- analyses_scheduled * nPatients
  analyses_nPatients
  efficacy_thresholds

  # Simulation Parameters
  nSims <- nSims
  pb <- txtProgressBar(
    min = 0, max = nSims,
    initial = 0, style = 3)
  trialnum <- numeric(nSims)

  or <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  lcl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  ucl <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  pvalue <- data.frame(matrix(ncol = nLooks, nrow = nSims))
  success <- data.frame(matrix(ncol = nLooks, nrow = nSims))

  # provide column names

  strings <- c("or_%d", "lcl_%d", "ucl_%d", "pvalue_%d", "success_%d")

  colnames(or) <- sprintf("or_%d", (1:nLooks))
  colnames(lcl) <- sprintf("lcl_%d", (1:nLooks))
  colnames(ucl) <- sprintf("ucl_%d", (1:nLooks))
  colnames(pvalue) <- sprintf("pvalue_%d", (1:nLooks))
  colnames(success) <- sprintf("success_%d", (1:nLooks))
  df <- data.frame(trialnum, or, lcl, ucl, pvalue, success, overall_success)
  overall_success <- numeric(nSims)


  set.seed(seed) # this sets the random seed for your results to be reproducible

  for (i in 1:nSims) {
    trialnum[i] <- i
    pid <- seq(1, nPatients, by = 1)
    treatment <- rep(1:2, nPatients / 2)
    deathprob <- numeric(nPatients)
    deathprob[treatment == 1] <- death1
    deathprob[treatment == 2] <- death2
    death <- rbinom(nPatients, 1, deathprob)
    trialdata <- data.frame(cbind(pid, treatment, death))

    for (j in 1:nLooks) {
      analysisdata <- subset(trialdata, pid <= analyses_nPatients[j])
      model <- glm(death ~ treatment, family = binomial(link = "logit"), data = analysisdata)
      or[i, j] <- exp(summary(model)$coefficients[2]) # compute odds ratio
      lcl[i, j] <- exp(confint.default((model))[2, 1]) # compute lower confidence limit
      ucl[i, j] <- exp(confint.default((model))[2, 2]) # compute upper confidence limit
      pvalue[i, j] <- summary(model)$coefficients[8]
      success[i, j] <- ifelse(or[i, j] < 1 & pvalue[i, j] < efficacy_thresholds[j], 1, 0)
    }

    overall_success[i] <- 0

    for (j in 1:nLooks)
    {
      if (success[i, j] == 1) {
        overall_success[i] <- 1
      }
    }
    setTxtProgressBar(pb, i)
  }

  df <- data.frame(trialnum, or, lcl, ucl, pvalue, success, overall_success)

  simulation_results <- data.frame(matrix(vector(), nrow = nPatients, ncol = (length(df))))
  colnames(simulation_results) <- c("trialnum", (do.call(rbind, lapply(1:length(strings),
    FUN = function(j) {
      (do.call(rbind, lapply(1:nLooks,
        FUN = function(i) ((sprintf(strings, i)))
      )[]))[, j]
    }
  ))), "overall_success")
  simulation_results[intersect(names(df), names(simulation_results))] <- df[intersect(names(df), names(simulation_results))]

  results <- list(
    summary(design),
    plot(design, 1), plot(designSet, type = 1),
    head(simulation_results, n = 10), table(overall_success),
    simulation_results
  )
  return(results)
}
