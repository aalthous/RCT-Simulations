# RCT-Simulation-v1

In this repository, I will provide base R code for simulations of randomized controlled trials.  

I will begin with fairly straightforward examples, though over time these posts may expand to more sophisticated code as time and interest allows.  For the first few examples I'm trying to minimize the need for advanced R knowledge (only installing packages if absolutely needed) in an effort to avoid hurdles for newer R users, and I'm annotating code as much as practical, but some basic knowledge of R (as well as general knowledge of RCT's and commonly performed statistical analyses) will facilitate your efforts to use the code stored here.

1. RCT_Binary_Outcome provides code to simulate two-group parallel-arm randomized trials with a binary outcome.
2. RCT_Binary_Outcome_Interim_Efficacy_Manual provides code to simulate two-group parallel-arm randomized trials with a binary outcome and one planned interim analysis with an efficacy stoppong rule.  The example code uses an O'Brien-Fleming approach with one interim analysis taking place at 50% observed information, using a threshold of p<0.0054 to stop at the interim analysis, otherwise using a threshold of p<0.0492 to declare success at the final analysis.
3. RCT_Binary_Outcome_Interim_Efficacy_rpact (work in progress)
4. RCT_Time_To_Event_Outcome provides code to simulate two-group parallel-arm randomized trials with a time-to-event outcome.
5. (planned) RCT_Time_To_Event_Outcome_Interim_Efficacy

More instructions will follow.
