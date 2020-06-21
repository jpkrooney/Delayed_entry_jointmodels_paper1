library(tidyverse)
library(JMbayes)
library(parallel)
library(ggfortify)
library(openxlsx)
#library(JMTools)

#####
# Define functions to help extract parameters and calculate AUC's
extract_sim_estimates <- function(JMlist, scenario, startM=1, model ){
    estimates <- list()
    for (i in 1:length(JMlist)){
        temp <- data.frame(Scenario = scenario, Model = model, M = startM + i -1,
                           rbind(summary(JMlist[[i]])$"CoefTable-Long",
                                 summary(JMlist[[i]])$"CoefTable-Event",
                                 D_1_1 = c( summary(JMlist[[i]])$D[1, 1], NA, NA, NA, NA, NA),
                                 D_1_2 = c( summary(JMlist[[i]])$D[1, 2], NA, NA, NA, NA, NA),
                                 D_2_1 = c( summary(JMlist[[i]])$D[2, 1], NA, NA, NA, NA, NA),
                                 D_2_2 = c( summary(JMlist[[i]])$D[2, 2], NA, NA, NA, NA, NA),
                                 DIC = c( summary(JMlist[[i]])$DIC, NA, NA, NA, NA, NA),
                                 pD = c( summary(JMlist[[i]])$pD, NA, NA, NA, NA, NA),
                                 LPML = c( summary(JMlist[[i]])$LPML, NA, NA, NA, NA, NA)
                           ))
        estimates[[i]] <- data.frame(Variable = row.names(temp), temp)
        
    }
    estimates_df <- do.call("rbind", estimates)
    row.names(estimates_df) <- NULL
    return(estimates_df)
}



# Load dataset 1
datasets1 <- readRDS("Data/sim_data1.RDS")

# Optional code to visualise survival curves from all data in a given scenario
# Plot raw survival curves of each dataset in dataset1 - can adapt to other datasets
#sfits <- list()
#sfits <- lapply(1:length(datasets1), function(x)
#                data.frame(M = x, 
#                fortify( survfit (Surv (surv_t, vital_st) ~ group, data = datasets1[[x]]$S_data) ) ) )
#sfits <- do.call("rbind", sfits)
#sfits$M <- factor(sfits$M)
#
#ggplot(sfits, aes(x = time, y = surv, group = M, col = strata)) +
#    theme(panel.grid=element_blank()) +
#    geom_point(size=0.1, alpha=0.5)
    

###### Extract estimates
# Scenario 1
JM_Scen1_modA <- readRDS("Results/Sim_Jms/JM_Scen1_modA.RDS")
scenario1_ests_a <- extract_sim_estimates(JM_Scen1_modA, scenario = 1, model = "A", startM = 1)
rm(JM_Scen1_modA)
gc()
JM_Scen1_modB <- readRDS("Results/Sim_Jms/JM_Scen1_modB.RDS")
scenario1_ests_b <- extract_sim_estimates(JM_Scen1_modB, scenario = 1, model = "B", startM = 1)
rm(JM_Scen1_modB)
gc()
JM_Scen1_modC <- readRDS("Results/Sim_Jms/JM_Scen1_modC.RDS")
scenario1_ests_c <- extract_sim_estimates(JM_Scen1_modC, scenario = 1, model = "C", startM = 1)
rm(JM_Scen1_modC)
gc()
JM_Scen1_modD <- readRDS("Results/Sim_Jms/JM_Scen1_modD.RDS")
scenario1_ests_d <- extract_sim_estimates(JM_Scen1_modD, scenario = 1, model = "D", startM = 1)
rm(JM_Scen1_modD)
gc()

# Scenario 2
JM_Scen2_modA <- readRDS("Results/Sim_Jms/JM_Scen2_modA.RDS")
scenario2_ests_a <- extract_sim_estimates(JM_Scen2_modA, scenario = 2, model = "A", startM = 1)
rm(JM_Scen2_modA)
gc()
JM_Scen2_modB <- readRDS("Results/Sim_Jms/JM_Scen2_modB.RDS")
scenario2_ests_b <- extract_sim_estimates(JM_Scen2_modB, scenario = 2, model = "B", startM = 1)
rm(JM_Scen2_modB)
gc()
JM_Scen2_modC <- readRDS("Results/Sim_Jms/JM_Scen2_modC.RDS")
scenario2_ests_c <- extract_sim_estimates(JM_Scen2_modC, scenario = 2, model = "C", startM = 1)
rm(JM_Scen2_modC)
gc()
JM_Scen2_modD <- readRDS("Results/Sim_Jms/JM_Scen2_modD.RDS")
scenario2_ests_d <- extract_sim_estimates(JM_Scen2_modD, scenario = 2, model = "D", startM = 1)
rm(JM_Scen2_modD)
gc()

# Scenario 3
JM_Scen3_modA <- readRDS("Results/Sim_Jms/JM_Scen3_modA.RDS")
scenario3_ests_a <- extract_sim_estimates(JM_Scen3_modA, scenario = 3, model = "A", startM = 1)
rm(JM_Scen3_modA)
gc()
JM_Scen3_modB <- readRDS("Results/Sim_Jms/JM_Scen3_modB.RDS")
scenario3_ests_b <- extract_sim_estimates(JM_Scen3_modB, scenario = 3, model = "B", startM = 1)
rm(JM_Scen3_modB)
gc()
JM_Scen3_modC <- readRDS("Results/Sim_Jms/JM_Scen3_modC.RDS")
scenario3_ests_c <- extract_sim_estimates(JM_Scen3_modC, scenario = 3, model = "C", startM = 1)
rm(JM_Scen3_modC)
gc()
JM_Scen3_modD <- readRDS("Results/Sim_Jms/JM_Scen3_modD.RDS")
scenario3_ests_d <- extract_sim_estimates(JM_Scen3_modD, scenario = 3, model = "D", startM = 1)
rm(JM_Scen3_modD)
gc()

# Scenario 4
JM_Scen4_modA <- readRDS("Results/Sim_Jms/JM_Scen4_modA.RDS")
scenario4_ests_a <- extract_sim_estimates(JM_Scen4_modA, scenario = 4, model = "A", startM = 1)
rm(JM_Scen4_modA)
gc()
JM_Scen4_modB <- readRDS("Results/Sim_Jms/JM_Scen4_modB.RDS")
scenario4_ests_b <- extract_sim_estimates(JM_Scen4_modB, scenario = 4, model = "B", startM = 1)
rm(JM_Scen4_modB)
gc()
JM_Scen4_modC <- readRDS("Results/Sim_Jms/JM_Scen4_modC.RDS")
scenario4_ests_c <- extract_sim_estimates(JM_Scen4_modC, scenario = 4, model = "C", startM = 1)
rm(JM_Scen4_modC)
gc()
JM_Scen4_modD <- readRDS("Results/Sim_Jms/JM_Scen4_modD.RDS")
scenario4_ests_d <- extract_sim_estimates(JM_Scen4_modD, scenario = 4, model = "D", startM = 1)
rm(JM_Scen4_modD)
gc()


# Scenario 5
JM_Scen5_modA <- readRDS("Results/Sim_Jms/JM_Scen5_modA.RDS")
scenario5_ests_a <- extract_sim_estimates(JM_Scen5_modA, scenario = 5, model = "A", startM = 1)
rm(JM_Scen5_modA)
gc()
JM_Scen5_modB <- readRDS("Results/Sim_Jms/JM_Scen5_modB.RDS")
scenario5_ests_b <- extract_sim_estimates(JM_Scen5_modB, scenario = 5, model = "B", startM = 1)
rm(JM_Scen5_modB)
gc()
JM_Scen5_modC <- readRDS("Results/Sim_Jms/JM_Scen5_modC.RDS")
scenario5_ests_c <- extract_sim_estimates(JM_Scen5_modC, scenario = 5, model = "C", startM = 1)
rm(JM_Scen5_modC)
gc()
JM_Scen5_modD <- readRDS("Results/Sim_Jms/JM_Scen5_modD.RDS")
scenario5_ests_d <- extract_sim_estimates(JM_Scen5_modD, scenario = 5, model = "D", startM = 1)
rm(JM_Scen5_modD)
gc()



# Combine estimates 

scenario1_ests <- rbind(scenario1_ests_a, scenario1_ests_b, scenario1_ests_c, scenario1_ests_d)
scenario2_ests <- rbind(scenario2_ests_a, scenario2_ests_b, scenario2_ests_c, scenario2_ests_d)
scenario3_ests <- rbind(scenario3_ests_a, scenario3_ests_b, scenario3_ests_c, scenario3_ests_d)
scenario4_ests <- rbind(scenario4_ests_a, scenario4_ests_b, scenario4_ests_c, scenario4_ests_d)
scenario5_ests <- rbind(scenario5_ests_a, scenario5_ests_b, scenario5_ests_c, scenario5_ests_d)


#
df1 <- rbind(scenario1_ests, scenario2_ests, scenario3_ests, scenario4_ests,
             scenario5_ests)

# Write all results to a file
write.csv(df1, "Results/Sim_results/All_model_estimates.csv", row.names = FALSE)
df1 <- read.csv("Results/Sim_results/All_model_estimates.csv")


#df1$Variable <- as.character(df1$Variable)

df1 <- df1 %>%
    mutate(Variable = case_when(Variable == "(Intercept)" ~ "Intercept",
                                Variable %in% c("time", "poly(time, 2)1") ~ "time",
                                Variable == "poly(time, 2)2" ~ "time^2",          
                                Variable %in% c("adj_time", "poly(adj_time, 2)1") ~ "adj_time",
                                Variable == "poly(adj_time, 2)2" ~ "adj_time^2",          
                                Variable == "entry_t2" ~ "Delayed entry Hazard",
                                Variable == "Assoct" ~ "Alpha",
                                Variable == "group2" ~ "Sub-group Rel. Hazard",
                                Variable == "D_1_1" ~ "Sigma_1_1",
                                Variable == "D_1_2" ~ "Sigma_1_2",
                                Variable == "D_2_1" ~ "Sigma_2_1",
                                Variable == "D_2_2" ~ "Sigma_2_2",
                                Variable == "DIC" ~ "DIC",
                                Variable == "pD" ~ "pD",
                                Variable == "tauBs" ~ "tauBs",
                                Variable == "LPML" ~ "LPML"
    ))
df1$Variable <- as.factor(df1$Variable)
df1$Scenario <- as.factor(df1$Scenario)
df1$Model <- as.factor(df1$Model)
df1$X2.5. <- as.numeric(df1$X2.5.)
df1$X97.5. <- as.numeric(df1$X97.5.)



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
df1 <- df1 %>% 
    mutate(True_value = case_when(Variable == "Alpha" ~ -0.05,
                                  Variable == "Sub-group Rel. Hazard" ~ log(1.5),
                                  Variable == "Intercept" & Model %in% c("A", "B") ~ 60,
                                  Variable == "Intercept" & Model %in% c("C", "D") &
                                      Scenario == 1 ~ scen1_mean_int,
                                  Variable == "Intercept" & Model %in% c("C", "D") &
                                      Scenario == 2 ~ scen2_mean_int,
                                  Variable == "Intercept" & Model %in% c("C", "D") &
                                      Scenario == 3 ~ scen3_mean_int,
                                  Variable == "Intercept" & Model %in% c("C", "D") &
                                      Scenario == 4 ~ scen4_mean_int,
                                  Variable == "Intercept" & Model %in% c("C", "D") &
                                      Scenario == 5 ~ scen5_mean_int,
                                  Variable == "time" & Scenario %in% c(1, 2, 4)~ -1,
                                  Variable == "adj_time" & Scenario %in% c(1, 2, 4)~ -1,
                                  Variable == "time" & Scenario ==3 ~ scen3_mean_slope,
                                  Variable == "adj_time" & Scenario ==3 ~ scen3_mean_slope,
                                  Variable == "Delayed entry Hazard" & Scenario == 1 ~ 0,
                                  Variable == "Delayed entry Hazard" & Scenario %in% c(2:5) ~ -0.05
    ))

# Tag each row as estimate within 95%CI or not
df1$cov <- ifelse(df1$X2.5. < df1$True_value & df1$X97.5. > df1$True_value, 1, 0)

# Calculate squared error between each estimate and truth for later MSE calculation
df1$Value <- as.numeric(df1$Value)
df1$sqr_err <- (df1$True_value - df1$Value) ^ 2



# Calculate mean estimates from simulations
tab_estimates <- df1[df1$Variable %in% c("Alpha", "Sub-group Rel. Hazard", "Intercept", "time",
                        "adj_time", "time^2", "adj_time^2", "Delayed entry Hazard",
                        "RE_Intercept", "RE_slope"), ] %>%
    group_by(Variable, Scenario, Model) %>%
    dplyr::summarize(Mean = mean(Value),
                     Coverage = mean(cov),
                     MSE  = mean(sqr_err),
                     True_value = mean(True_value))


# Calculate bias
tab_estimates <- tab_estimates %>% 
    mutate(Bias = Mean - True_value)

# Calculate bias%
tab_estimates <- tab_estimates %>% 
    mutate(Bias_pct = Bias * 100 / True_value)


# Reorder columns
tab_estimates <- dplyr::select(tab_estimates, Scenario, Variable, Model, True_value,
                        Mean, MSE, Bias, Bias_pct, Coverage ) %>% 
                        arrange(Scenario, Variable, Model)

# Write results to a file
write.csv(tab_estimates, "Results/Sim_results/Estimates_summaries_new.csv", row.names = FALSE)

# Write results to a worksheet for each scenario
## Create a new workbook
wb1 <- createWorkbook()

for(i in 1:5){
    # Add worksheet
    name <- paste0("Scenario", i)
    addWorksheet(wb1, name)
    writeData(wb1, sheet = name, tab_estimates %>% filter(Scenario == i))
}
# Save workbook to file
saveWorkbook(wb1, "Results/Sim_results/JM_ests_by_scenario.xlsx", overwrite = TRUE)

