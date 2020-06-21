library(tidyverse)
library(broom)
library(parallel)
library(ggfortify)


#####
# Define functions to extract parameters 
extract_simCox_estimates <- function(coxlist, scenario, startM=1, model ){
    estimates <- list()
    for (i in 1:length(coxlist)){
        estimates[[i]] <- data.frame(Scenario = scenario, Model = model, M = startM + i -1,
                           rbind(tidy(coxlist[[i]], exponentiate = FALSE, )))

    }
    estimates_df <- do.call("rbind", estimates)
    return(estimates_df)
}



# Load Cox models Scenario 1
coxScen1_modA <- readRDS("Results/Sim_sub-models/coxScen1_modA.RDS")
coxScen1_modB <- readRDS("Results/Sim_sub-models/coxScen1_modB.RDS")
coxScen1_modC <- readRDS("Results/Sim_sub-models/coxScen1_modC.RDS")
coxScen1_modD <- readRDS("Results/Sim_sub-models/coxScen1_modD.RDS")
# Extract estimates
ests_cox_scen1ModA <- extract_simCox_estimates(coxScen1_modA, scenario = 1, model = "A")
ests_cox_scen1ModB <- extract_simCox_estimates(coxScen1_modB, scenario = 1, model = "B")
ests_cox_scen1ModC <- extract_simCox_estimates(coxScen1_modC, scenario = 1, model = "C")
ests_cox_scen1ModD <- extract_simCox_estimates(coxScen1_modD, scenario = 1, model = "D")

# Load Cox models Scenario 2
coxScen2_modA <- readRDS("Results/Sim_sub-models/coxScen2_modA.RDS")
coxScen2_modB <- readRDS("Results/Sim_sub-models/coxScen2_modB.RDS")
coxScen2_modC <- readRDS("Results/Sim_sub-models/coxScen2_modC.RDS")
coxScen2_modD <- readRDS("Results/Sim_sub-models/coxScen2_modD.RDS")
# Extract estimates
ests_cox_scen2ModA <- extract_simCox_estimates(coxScen2_modA, scenario = 2, model = "A")
ests_cox_scen2ModB <- extract_simCox_estimates(coxScen2_modB, scenario = 2, model = "B")
ests_cox_scen2ModC <- extract_simCox_estimates(coxScen2_modC, scenario = 2, model = "C")
ests_cox_scen2ModD <- extract_simCox_estimates(coxScen2_modD, scenario = 2, model = "D")

# Load Cox models Scenario 3
coxScen3_modA <- readRDS("Results/Sim_sub-models/coxScen3_modA.RDS")
coxScen3_modB <- readRDS("Results/Sim_sub-models/coxScen3_modB.RDS")
coxScen3_modC <- readRDS("Results/Sim_sub-models/coxScen3_modC.RDS")
coxScen3_modD <- readRDS("Results/Sim_sub-models/coxScen3_modD.RDS")
# Extract estimates
ests_cox_scen3ModA <- extract_simCox_estimates(coxScen3_modA, scenario = 3, model = "A")
ests_cox_scen3ModB <- extract_simCox_estimates(coxScen3_modB, scenario = 3, model = "B")
ests_cox_scen3ModC <- extract_simCox_estimates(coxScen3_modC, scenario = 3, model = "C")
ests_cox_scen3ModD <- extract_simCox_estimates(coxScen3_modD, scenario = 3, model = "D")

# Load Cox models Scenario 4
coxScen4_modA <- readRDS("Results/Sim_sub-models/coxScen4_modA.RDS")
coxScen4_modB <- readRDS("Results/Sim_sub-models/coxScen4_modB.RDS")
coxScen4_modC <- readRDS("Results/Sim_sub-models/coxScen4_modC.RDS")
coxScen4_modD <- readRDS("Results/Sim_sub-models/coxScen4_modD.RDS")
# Extract estimates
ests_cox_scen4ModA <- extract_simCox_estimates(coxScen4_modA, scenario = 4, model = "A")
ests_cox_scen4ModB <- extract_simCox_estimates(coxScen4_modB, scenario = 4, model = "B")
ests_cox_scen4ModC <- extract_simCox_estimates(coxScen4_modC, scenario = 4, model = "C")
ests_cox_scen4ModD <- extract_simCox_estimates(coxScen4_modD, scenario = 4, model = "D")

# Load Cox models Scenario 5
coxScen5_modA <- readRDS("Results/Sim_sub-models/coxScen5_modA.RDS")
coxScen5_modB <- readRDS("Results/Sim_sub-models/coxScen5_modB.RDS")
coxScen5_modC <- readRDS("Results/Sim_sub-models/coxScen5_modC.RDS")
coxScen5_modD <- readRDS("Results/Sim_sub-models/coxScen5_modD.RDS")
# Extract estimates
ests_cox_scen5ModA <- extract_simCox_estimates(coxScen5_modA, scenario = 5, model = "A")
ests_cox_scen5ModB <- extract_simCox_estimates(coxScen5_modB, scenario = 5, model = "B")
ests_cox_scen5ModC <- extract_simCox_estimates(coxScen5_modC, scenario = 5, model = "C")
ests_cox_scen5ModD <- extract_simCox_estimates(coxScen5_modD, scenario = 5, model = "D")


ests_cox_All <- rbind(ests_cox_scen1ModA, ests_cox_scen1ModB, ests_cox_scen1ModC, ests_cox_scen1ModD,
                      ests_cox_scen2ModA, ests_cox_scen2ModB, ests_cox_scen2ModC, ests_cox_scen2ModD,
                      ests_cox_scen3ModA, ests_cox_scen3ModB, ests_cox_scen3ModC, ests_cox_scen3ModD,
                      ests_cox_scen4ModA, ests_cox_scen4ModB, ests_cox_scen4ModC, ests_cox_scen4ModD,
                      ests_cox_scen5ModA, ests_cox_scen5ModB, ests_cox_scen5ModC, ests_cox_scen5ModD
                      )

ests_cox_All$term <- ifelse( ests_cox_All$term == "group2", "Sub-group Rel. Hazard",
                             ifelse(ests_cox_All$term == "entry_t2", "Delayed entry Hazard", NA))


# Assign true values to df1 for later coverage and MSE calculations
dfres <- ests_cox_All %>% 
    mutate(True_value = case_when(term == "Sub-group Rel. Hazard" ~ log(1.5),
                                  term == "Delayed entry Hazard" & Scenario == 1 ~ 0,
                                  term == "Delayed entry Hazard" & Scenario %in% c(2:5) ~ -0.05
    ))

# Tag each row as estimate within 95%CI or not
dfres$cover <- ifelse(dfres$conf.low < dfres$True_value & dfres$conf.high > dfres$True_value, 1, 0)

# Calculate squared error between each estimate and truth for later MSE calculation
dfres$sqr_err <- (dfres$True_value - dfres$estimate) ^ 2


# Calculate mean estimates from simulations
tab_estimates <- dfres %>%
    group_by(term, Scenario, Model) %>%
    dplyr::summarize(Mean = mean(estimate),
                     Coverage = mean(cover),
                     MSE  = mean(sqr_err),
                     True_value = mean(True_value))  # Each value of True value shoud be the same with a defined group

# Calculate bias
tab_estimates <- tab_estimates %>% 
    mutate(Bias = Mean - True_value)

# Calculate bias%
tab_estimates <- tab_estimates %>% 
    mutate(Bias_pct = Bias * 100 / True_value)

# Reorder columns and rows
tab_estimates <- dplyr::select(tab_estimates, term, Scenario, Model, True_value,
                        Mean, MSE, Bias, Bias_pct, Coverage ) %>% 
                        arrange(Scenario, term, Model)


# Write results to a file
write.csv(tab_estimates, "Results/Sim_results/Cox_ests_summaries.csv", row.names = FALSE)



##### Linear MM estimates

# Define functions to extract LMM fixed parameters
extract_LMM_fixef <- function(LMMlist, scenario, startM=1, model ){
    estimates <- list()
    for (i in 1:length(LMMlist)){
        estimates[[i]] <- data.frame(Scenario = scenario, Model = model, M = startM + i -1,
                                     rbind(tidy(LMMlist[[i]], effects = "fixed" )))
        
    }
    estimates_df <- do.call("rbind", estimates)
    return(estimates_df)
}






# Load Cox models Scenario 1
LmmScen1_time <- readRDS("Results/Sim_sub-models/lmeScen1_time.RDS")
LmmScen1_adjtime <- readRDS("Results/Sim_sub-models/lmeScen1_adjtime.RDS")

# Load Cox models Scenario 2
LmmScen2_time <- readRDS("Results/Sim_sub-models/lmeScen2_time.RDS")
LmmScen2_adjtime <- readRDS("Results/Sim_sub-models/lmeScen2_adjtime.RDS")

# Load Cox models Scenario 3
LmmScen3_time <- readRDS("Results/Sim_sub-models/lmeScen3_time.RDS")
LmmScen3_adjtime <- readRDS("Results/Sim_sub-models/lmeScen3_adjtime.RDS")

# Load Cox models Scenario 4
LmmScen4_time <- readRDS("Results/Sim_sub-models/lmeScen4_time.RDS")
LmmScen4_adjtime <- readRDS("Results/Sim_sub-models/lmeScen4_adjtime.RDS")

# Load Cox models Scenario 5
LmmScen5_time <- readRDS("Results/Sim_sub-models/lmeScen5_time.RDS")
LmmScen5_adjtime <- readRDS("Results/Sim_sub-models/lmeScen5_adjtime.RDS")


# Extract estimates
ests_LMM_scen1_time <- extract_LMM_fixef(LmmScen1_time, scenario = 1, model = "time")
ests_LMM_scen2_time <- extract_LMM_fixef(LmmScen2_time, scenario = 2, model = "time")
ests_LMM_scen3_time <- extract_LMM_fixef(LmmScen3_time, scenario = 3, model = "time")
ests_LMM_scen4_time <- extract_LMM_fixef(LmmScen4_time, scenario = 4, model = "time")
ests_LMM_scen5_time <- extract_LMM_fixef(LmmScen5_time, scenario = 5, model = "time")

ests_LMM_scen1_adjtime <- extract_LMM_fixef(LmmScen1_adjtime, scenario = 1, model = "adjtime")
ests_LMM_scen2_adjtime <- extract_LMM_fixef(LmmScen2_adjtime, scenario = 2, model = "adjtime")
ests_LMM_scen3_adjtime <- extract_LMM_fixef(LmmScen3_adjtime, scenario = 3, model = "adjtime")
ests_LMM_scen4_adjtime <- extract_LMM_fixef(LmmScen4_adjtime, scenario = 4, model = "adjtime")
ests_LMM_scen5_adjtime <- extract_LMM_fixef(LmmScen5_adjtime, scenario = 5, model = "adjtime")

# 
ests_LMM_All <- rbind(ests_LMM_scen1_time, ests_LMM_scen2_time, ests_LMM_scen3_time,
                      ests_LMM_scen4_time, ests_LMM_scen5_time,
                      ests_LMM_scen1_adjtime, ests_LMM_scen2_adjtime, ests_LMM_scen3_adjtime,
                      ests_LMM_scen4_adjtime, ests_LMM_scen5_adjtime)

ests_LMM_All$term <- ifelse( ests_LMM_All$term == "(Intercept)", "Intercept",
                         ifelse( ests_LMM_All$term %in% c("time", "adj_time"),
                                 "Slope", ests_LMM_All$term))


# For all scenarios, true values for the intercepts on the adjusted timeline (i.e. t_adj = 0 at time of diagnosis) can be calculated from the individualised b0, b1 and delayed entry times.

datasets1 <- readRDS("Data/sim_data1.RDS")
obs_ints_scen1 <- unlist( lapply(1: length(datasets1), function(x)
    datasets1[[x]]$L_data$B0 + ( datasets1[[x]]$L_data$B1 * datasets1[[x]]$L_data$entry_t  ) ) )
trunc_ints_scen1 <- unlist( lapply(1: length(datasets1), function(x)
    datasets1[[x]]$trunc_Ldata$B0 + ( datasets1[[x]]$trunc_Ldata$B1 * datasets1[[x]]$trunc_Ldata$entry_t
    )))
scen1_mean_int <- mean( c(obs_ints_scen1, trunc_ints_scen1) )
rm(datasets1)
#
datasets2 <- readRDS("Data/sim_data2.RDS")
obs_ints_scen2 <- unlist( lapply(1: length(datasets2), function(x)
    datasets2[[x]]$L_data$B0 + ( datasets2[[x]]$L_data$B1 * datasets2[[x]]$L_data$entry_t  ) ) )
trunc_ints_scen2 <- unlist( lapply(1: length(datasets2), function(x)
    datasets2[[x]]$trunc_Ldata$B0 + ( datasets2[[x]]$trunc_Ldata$B1 * datasets2[[x]]$trunc_Ldata$entry_t
    )))
scen2_mean_int <- mean( c(obs_ints_scen2, trunc_ints_scen2) )
rm(datasets2)

# Note for dataset 3 must also calculate mean slope
datasets3 <- readRDS("Data/sim_data3.RDS")
obs_ints_scen3 <- unlist( lapply(1: length(datasets3), function(x)
    datasets3[[x]]$L_data$B0 + ( datasets3[[x]]$L_data$B1 * datasets3[[x]]$L_data$entry_t  ) ) )
trunc_ints_scen3 <- unlist( lapply(1: length(datasets3), function(x)
    datasets3[[x]]$trunc_Ldata$B0 + ( datasets3[[x]]$trunc_Ldata$B1 * datasets3[[x]]$trunc_Ldata$entry_t
    )))
scen3_mean_int <- mean( c(obs_ints_scen3, trunc_ints_scen3) )
#
obs_slope_scen3 <- unlist( lapply(1: length(datasets3), function(x) datasets3[[x]]$S_data$B1 ) )
trunc_slope_scen3 <- unlist( lapply(1: length(datasets3), function(x) datasets3[[x]]$trunc_Sdata$B1 ) )
scen3_mean_slope <- mean( c(obs_slope_scen3, trunc_slope_scen3 ) )

# Note for dataset 4 there is also an effect of entry_t on Y ltrunc_long_beta = 0.1 * dx_delay
datasets4 <- readRDS("Data/sim_data4.RDS")

obs_ints_scen4 <- unlist( lapply(1: length(datasets4), function(x)
    datasets4[[x]]$L_data$B0 + ( datasets4[[x]]$L_data$B1 * datasets4[[x]]$L_data$entry_t  ) +
        0.1 * datasets4[[x]]$L_data$entry_t) )
trunc_ints_scen4 <- unlist( lapply(1: length(datasets4), function(x)
    datasets4[[x]]$trunc_Ldata$B0 + ( datasets4[[x]]$trunc_Ldata$B1 * datasets4[[x]]$trunc_Ldata$entry_t  ) +
        0.1 * datasets4[[x]]$trunc_Ldata$entry_t) )
scen4_mean_int <- mean( c(obs_ints_scen4, trunc_ints_scen4) )
rm(datasets4)
#
datasets5 <- readRDS("Data/sim_data5.RDS")
obs_ints_scen5 <- unlist( lapply(1: length(datasets5), function(x)
    datasets5[[x]]$L_data$B0 + ( datasets5[[x]]$L_data$B1 * datasets5[[x]]$L_data$entry_t  ) ) )
trunc_ints_scen5 <- unlist( lapply(1: length(datasets5), function(x)
    datasets5[[x]]$trunc_Ldata$B0 + ( datasets5[[x]]$trunc_Ldata$B1 * datasets5[[x]]$trunc_Ldata$entry_t
    )))
scen5_mean_int <- mean( c(obs_ints_scen5, trunc_ints_scen5) )
rm(datasets5)




# Assign true values to df1 for later coverage and MSE calculations
#dfresLMM <- ests_LMM_All %>% 
#    mutate(True_value = case_when(term == "Intercept" ~ 60,
#                                  term =="Slope" & Scenario %in% c(1, 2, 4) ~ -1,
#                                  term =="Slope" & Scenario %in% c(3) ~ scen3_mean_slope
#    ))
# Assign true values to df1 for later coverage and MSE calculations
dfresLMM <- ests_LMM_All %>% 
    mutate(True_value = case_when(term == "Intercept" & Model == "time" ~ 60,
                                  term == "Intercept" & Model == "adjtime" &
                                      Scenario == 1 ~ scen1_mean_int,
                                  term == "Intercept" & Model =="adjtime" &
                                      Scenario == 2 ~ scen2_mean_int,
                                  term == "Intercept" & Model == "adjtime" &
                                      Scenario == 3 ~ scen3_mean_int,
                                  term == "Intercept" & Model == "adjtime" &
                                      Scenario == 4 ~ scen4_mean_int,
                                  term == "Intercept" & Model == "adjtime" &
                                      Scenario == 5 ~ scen5_mean_int,
                                  term == "Slope" & Scenario %in% c(1, 2, 4)~ -1,
                                  term == "Slope" & Scenario ==3 ~ scen3_mean_slope))


# Calculate squared error between each estimate and truth for later MSE calculation
dfresLMM $sqr_err <- (dfresLMM$True_value - dfresLMM$estimate) ^ 2


tab_LMMests <- dfresLMM %>%
    group_by(term, Scenario, Model) %>%
    dplyr::summarize(Mean = mean(estimate),
                     MSE  = mean(sqr_err),
                     True_value = mean(True_value))  # Each value of True value shoud be the same with a defined group


# Calculate bias
tab_LMMests <- tab_LMMests %>% 
    mutate(Bias = Mean - True_value)

# Calculate bias%
tab_LMMests <- tab_LMMests %>% 
    mutate(Bias_pct = Bias * 100 / True_value)

# Reorder columns and rows
tab_LMMests <- dplyr::select(tab_LMMests, Scenario, term, Model, True_value,
                        Mean, MSE, Bias, Bias_pct) %>% 
    arrange(Scenario, Model, term)


# Write results to a file
write.csv(tab_LMMests, "Results/Sim_results/LMM_ests_summaries.csv", row.names = FALSE)



