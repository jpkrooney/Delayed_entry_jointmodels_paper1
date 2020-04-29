library(tidyverse)
library(JMbayes)
library(future.apply)

#library(JMTools)

# Set plan mode and num workers for future.apply parallel processing
# Note each worker needs up to 4GB of RAM
plan(multisession, workers = 12, gc=TRUE)

# Set memory limit for futures to 1GB
options(future.globals.maxSize = 1000*1024^2)

# Set RNG algorithm to support parallel reproducibility
RNGkind(kind = "L'Ecuyer-CMRG")


# Load simulated datasets 1
datasets1 <- readRDS("Data/sim_data1.RDS")


# Read number of models from data
M <- length(datasets1)


###### Prepare cox and linear mixed models for scenario 1
# Make lists to hold submodel results
coxScen1_modA <- list(); coxScen1_modB <- list(); coxScen1_modC <- list(); coxScen1_modD <- list()
lmeScen1_time <- list(); lmeScen1_adjtime <- list()


# set a seed
set.seed(45623234)
for (i in 1: M){
    # Set timefix = FALSE to avoid coxph floating point error than can occur with ties
    # in simulated datasets as per
    # https://cran.r-project.org/web/packages/survival/vignettes/tiedtimes.pdf
    coxScen1_modA[[i]] <- coxph (Surv (surv_t, vital_st) ~ group + entry_t2,
                            data = datasets1[[i]]$S_data, x = TRUE, model = TRUE,
                            timefix = FALSE)
    coxScen1_modB[[i]] <- coxph (Surv (entry_t, surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets1[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen1_modC[[i]] <- coxph (Surv (adj_surv_t, vital_st) ~ group,
                                 data = datasets1[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen1_modD[[i]] <- coxph (Surv (adj_surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets1[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    
    # Longitudinal model on natural and adjusted time
    lmeScen1_time[[i]] <- lme(Y ~ time,
                         random = ~ time | ID, data = datasets1[[i]]$L_data,
                         control = lmeControl(opt = "optim"))
    lmeScen1_adjtime[[i]] <- lme(Y ~ adj_time,
                         random = ~ adj_time | ID, data = datasets1[[i]]$L_data,
                         control = lmeControl(opt = "optim")) 
}

### Save Cox and LME models to file
saveRDS(coxScen1_modA, "Results/Sim_sub-models/coxScen1_modA.RDS")
saveRDS(coxScen1_modB, "Results/Sim_sub-models/coxScen1_modB.RDS")
saveRDS(coxScen1_modC, "Results/Sim_sub-models/coxScen1_modC.RDS")
saveRDS(coxScen1_modD, "Results/Sim_sub-models/coxScen1_modD.RDS")
saveRDS(lmeScen1_time, "Results/Sim_sub-models/lmeScen1_time.RDS")
saveRDS(lmeScen1_adjtime, "Results/Sim_sub-models/lmeScen1_adjtime.RDS")


# Note seed will be specified for ach run of future_lapply which manages the
# L'Ecuyer-CMRGbetween algorithm internally whne provided with a seed 

# Build joint models for scenario 1 model A all datasets
JM_Scen1_modA <- future_lapply(1:M, function(x)
                        jointModelBayes(lmeScen1_time[[x]], coxScen1_modA[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
                    future.seed = 23578726)
# Save this batch
saveRDS(JM_Scen1_modA, "Results/Sim_Jms/JM_Scen1_modA.RDS")
# Clean-up to free up RAM
rm(JM_Scen1_modA)
gc()

###

# Build joint models for scenario 1 model B all datasets
JM_Scen1_modB <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen1_time[[x]], coxScen1_modB[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 3549734)
# Save this batch
saveRDS(JM_Scen1_modB, "Results/Sim_Jms/JM_Scen1_modB.RDS")
# Clean-up to free up RAM
rm(JM_Scen1_modB)
gc()

###

# Build joint models for scenario 1 model C all datasets
JM_Scen1_modC <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen1_adjtime[[x]], coxScen1_modC[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 60947309)
# Save this batch
saveRDS(JM_Scen1_modC, "Results/Sim_Jms/JM_Scen1_modC.RDS")
# Clean-up to free up RAM
rm(JM_Scen1_modC)
gc()

###

# Build joint models for scenario 1 model D all datasets
JM_Scen1_modD <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen1_adjtime[[x]], coxScen1_modD[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 4353950)
# Save this batch
saveRDS(JM_Scen1_modD, "Results/Sim_Jms/JM_Scen1_modD.RDS")
# Clean-up to free up RAM
rm(JM_Scen1_modD)
rm(datasets1)
gc()



##### Scenario 2

# Load simulated datasets 2
datasets2 <- readRDS("Data/sim_data2.RDS")

# Read number of models from data
M <- length(datasets2)


###### Prepare cox and linear mixed models for scenario 2
# Make lists to hold submodel results
coxScen2_modA <- list(); coxScen2_modB <- list(); coxScen2_modC <- list(); coxScen2_modD <- list()
lmeScen2_time <- list(); lmeScen2_adjtime <- list()


# set a seed
set.seed(3957659)
for (i in 1: M){
    # Set timefix = FALSE to avoid coxph floating point error than can occur with ties
    # in simulated datasets as per
    # https://cran.r-project.org/web/packages/survival/vignettes/tiedtimes.pdf
    coxScen2_modA[[i]] <- coxph (Surv (surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets2[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen2_modB[[i]] <- coxph (Surv (entry_t, surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets2[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen2_modC[[i]] <- coxph (Surv (adj_surv_t+.001, vital_st) ~ group,
                                 data = datasets2[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen2_modD[[i]] <- coxph (Surv (adj_surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets2[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    
    # Longitudinal model on natural and adjusted time
    lmeScen2_time[[i]] <- lme(Y ~ time,
                              random = ~ time | ID, data = datasets2[[i]]$L_data,
                              control = lmeControl(opt = "optim"))
    lmeScen2_adjtime[[i]] <- lme(Y ~ adj_time,
                                 random = ~ adj_time | ID, data = datasets2[[i]]$L_data,
                                 control = lmeControl(opt = "optim")) 
}

### Save Cox and LME models to file
saveRDS(coxScen2_modA, "Results/Sim_sub-models/coxScen2_modA.RDS")
saveRDS(coxScen2_modB, "Results/Sim_sub-models/coxScen2_modB.RDS")
saveRDS(coxScen2_modC, "Results/Sim_sub-models/coxScen2_modC.RDS")
saveRDS(coxScen2_modD, "Results/Sim_sub-models/coxScen2_modD.RDS")
saveRDS(lmeScen2_time, "Results/Sim_sub-models/lmeScen2_time.RDS")
saveRDS(lmeScen2_adjtime, "Results/Sim_sub-models/lmeScen2_adjtime.RDS")


# Note seed will be specified for ach run of future_lapply which manages the
# L'Ecuyer-CMRGbetween algorithm internally whne provided with a seed 

# Build joint models for scenario 2 model A all datasets
JM_Scen2_modA <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen2_time[[x]], coxScen2_modA[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 47596)
# Save this batch
saveRDS(JM_Scen2_modA, "Results/Sim_Jms/JM_Scen2_modA.RDS")
# Clean-up to free up RAM
rm(JM_Scen2_modA)
gc()

###

# Build joint models for scenario 2 model B all datasets
JM_Scen2_modB <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen2_time[[x]], coxScen2_modB[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 5906720)
# Save this batch
saveRDS(JM_Scen2_modB, "Results/Sim_Jms/JM_Scen2_modB.RDS")
# Clean-up to free up RAM
rm(JM_Scen2_modB)
gc()

###

# Build joint models for scenario 2 model C all datasets
JM_Scen2_modC <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen2_adjtime[[x]], coxScen2_modC[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 234586)
# Save this batch
saveRDS(JM_Scen2_modC, "Results/Sim_Jms/JM_Scen2_modC.RDS")
# Clean-up to free up RAM
rm(JM_Scen2_modC)
gc()

###

# Build joint models for scenario 2 model D all datasets
JM_Scen2_modD <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen2_adjtime[[x]], coxScen2_modD[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 2093850)
# Save this batch
saveRDS(JM_Scen2_modD, "Results/Sim_Jms/JM_Scen2_modD.RDS")
# Clean-up to free up RAM
rm(JM_Scen2_modD)
gc()



##### Scenario 3

# Load simulated datasets 3
datasets3 <- readRDS("Data/sim_data3.RDS")

# Read number of models from data
M <- length(datasets3)

###### Prepare cox and linear mixed models for scenario 3
# Make lists to hold submodel results
coxScen3_modA <- list(); coxScen3_modB <- list(); coxScen3_modC <- list(); coxScen3_modD <- list()
lmeScen3_time <- list(); lmeScen3_adjtime <- list()


# set a seed
set.seed(7985873)
for (i in 1: M){
    # Set timefix = FALSE to avoid coxph floating point error than can occur with ties
    # in simulated datasets as per
    # https://cran.r-project.org/web/packages/survival/vignettes/tiedtimes.pdf
    coxScen3_modA[[i]] <- coxph (Surv (surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets3[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen3_modB[[i]] <- coxph (Surv (entry_t, surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets3[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen3_modC[[i]] <- coxph (Surv (adj_surv_t+.001, vital_st) ~ group,
                                 data = datasets3[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen3_modD[[i]] <- coxph (Surv (adj_surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets3[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    
    # Longitudinal model on natural and adjusted time
    lmeScen3_time[[i]] <- lme(Y ~ time,
                              random = ~ time | ID, data = datasets3[[i]]$L_data,
                              control = lmeControl(opt = "optim"))
    lmeScen3_adjtime[[i]] <- lme(Y ~ adj_time,
                                 random = ~ adj_time | ID, data = datasets3[[i]]$L_data,
                                 control = lmeControl(opt = "optim")) 
}

### Save Cox and LME models to file
saveRDS(coxScen3_modA, "Results/Sim_sub-models/coxScen3_modA.RDS")
saveRDS(coxScen3_modB, "Results/Sim_sub-models/coxScen3_modB.RDS")
saveRDS(coxScen3_modC, "Results/Sim_sub-models/coxScen3_modC.RDS")
saveRDS(coxScen3_modD, "Results/Sim_sub-models/coxScen3_modD.RDS")
saveRDS(lmeScen3_time, "Results/Sim_sub-models/lmeScen3_time.RDS")
saveRDS(lmeScen3_adjtime, "Results/Sim_sub-models/lmeScen3_adjtime.RDS")


# Note seed will be specified for ach run of future_lapply which manages the
# L'Ecuyer-CMRGbetween algorithm internally whne provided with a seed 

# Build joint models for scenario 3 model A all datasets
JM_Scen3_modA <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen3_time[[x]], coxScen3_modA[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 33490560)
# Save this batch
saveRDS(JM_Scen3_modA, "Results/Sim_Jms/JM_Scen3_modA.RDS")
# Clean-up to free up RAM
rm(JM_Scen3_modA)
gc()

###

# Build joint models for scenario 3 model B all datasets
JM_Scen3_modB <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen3_time[[x]], coxScen3_modB[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 9284505)
# Save this batch
saveRDS(JM_Scen3_modB, "Results/Sim_Jms/JM_Scen3_modB.RDS")
# Clean-up to free up RAM
rm(JM_Scen3_modB)
gc()

###

# Build joint models for scenario 3 model C all datasets
JM_Scen3_modC <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen3_adjtime[[x]], coxScen3_modC[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 10239405)
# Save this batch
saveRDS(JM_Scen3_modC, "Results/Sim_Jms/JM_Scen3_modC.RDS")
# Clean-up to free up RAM
rm(JM_Scen3_modC)
gc()

###

# Build joint models for scenario 3 model D all datasets
JM_Scen3_modD <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen3_adjtime[[x]], coxScen3_modD[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 4553353)
# Save this batch
saveRDS(JM_Scen3_modD, "Results/Sim_Jms/JM_Scen3_modD.RDS")
# Clean-up to free up RAM
rm(JM_Scen3_modD)
gc()





##### Scenario 4

# Load simulated datasets 4
datasets4 <- readRDS("Data/sim_data4.RDS")

# Read number of models from data
M <- length(datasets4)

###### Prepare cox and linear mixed models for scenario 4
# Make lists to hold submodel results
coxScen4_modA <- list(); coxScen4_modB <- list(); coxScen4_modC <- list(); coxScen4_modD <- list()
lmeScen4_time <- list(); lmeScen4_adjtime <- list()


# set a seed
set.seed(6344783)
for (i in 1: M){
    # Set timefix = FALSE to avoid coxph floating point error than can occur with ties
    # in simulated datasets as per
    # https://cran.r-project.org/web/packages/survival/vignettes/tiedtimes.pdf
    coxScen4_modA[[i]] <- coxph (Surv (surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets4[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen4_modB[[i]] <- coxph (Surv (entry_t, surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets4[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen4_modC[[i]] <- coxph (Surv (adj_surv_t, vital_st) ~ group,
                                 data = datasets4[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen4_modD[[i]] <- coxph (Surv (adj_surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets4[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    
    # Longitudinal model on natural and adjusted time
    lmeScen4_time[[i]] <- lme(Y ~ time,
                              random = ~ time | ID, data = datasets4[[i]]$L_data,
                              control = lmeControl(opt = "optim"))
    lmeScen4_adjtime[[i]] <- lme(Y ~ adj_time,
                                 random = ~ adj_time | ID, data = datasets4[[i]]$L_data,
                                 control = lmeControl(opt = "optim")) 
}

### Save Cox and LME models to file
saveRDS(coxScen4_modA, "Results/Sim_sub-models/coxScen4_modA.RDS")
saveRDS(coxScen4_modB, "Results/Sim_sub-models/coxScen4_modB.RDS")
saveRDS(coxScen4_modC, "Results/Sim_sub-models/coxScen4_modC.RDS")
saveRDS(coxScen4_modD, "Results/Sim_sub-models/coxScen4_modD.RDS")
saveRDS(lmeScen4_time, "Results/Sim_sub-models/lmeScen4_time.RDS")
saveRDS(lmeScen4_adjtime, "Results/Sim_sub-models/lmeScen4_adjtime.RDS")


# Note seed will be specified for ach run of future_lapply which manages the
# L'Ecuyer-CMRGbetween algorithm internally whne provided with a seed 

# Build joint models for scenario 4 model A all datasets
JM_Scen4_modA <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen4_time[[x]], coxScen4_modA[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 5623345)
# Save this batch
saveRDS(JM_Scen4_modA, "Results/Sim_Jms/JM_Scen4_modA.RDS")
# Clean-up to free up RAM
rm(JM_Scen4_modA)
gc()

###

# Build joint models for scenario 4 model B all datasets
JM_Scen4_modB <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen4_time[[x]], coxScen4_modB[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 839562)
# Save this batch
saveRDS(JM_Scen4_modB, "Results/Sim_Jms/JM_Scen4_modB.RDS")
# Clean-up to free up RAM
rm(JM_Scen4_modB)
gc()

###

# Build joint models for scenario 4 model C all datasets
JM_Scen4_modC <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen4_adjtime[[x]], coxScen4_modC[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 34534)
# Save this batch
saveRDS(JM_Scen4_modC, "Results/Sim_Jms/JM_Scen4_modC.RDS")
# Clean-up to free up RAM
rm(JM_Scen4_modC)
gc()

###

# Build joint models for scenario 4 model D all datasets
JM_Scen4_modD <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen4_adjtime[[x]], coxScen4_modD[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 982956)
# Save this batch
saveRDS(JM_Scen4_modD, "Results/Sim_Jms/JM_Scen4_modD.RDS")
# Clean-up to free up RAM
rm(JM_Scen4_modD)
gc()





##### Scenario 5

# Load simulated datasets 5
datasets5 <- readRDS("Data/sim_data5.RDS")

# Read number of models from data
M <- length(datasets5)

###### Prepare cox and linear mixed models for scenario 5
# Make lists to hold submodel results
coxScen5_modA <- list(); coxScen5_modB <- list(); coxScen5_modC <- list(); coxScen5_modD <- list()
lmeScen5_time <- list(); lmeScen5_adjtime <- list()


# set a seed
set.seed(6355783)
for (i in 1: M){
    # Set timefix = FALSE to avoid coxph floating point error than can occur with ties
    # in simulated datasets as per
    # https://cran.r-project.org/web/packages/survival/vignettes/tiedtimes.pdf
    coxScen5_modA[[i]] <- coxph (Surv (surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets5[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen5_modB[[i]] <- coxph (Surv (entry_t, surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets5[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen5_modC[[i]] <- coxph (Surv (adj_surv_t, vital_st) ~ group,
                                 data = datasets5[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    coxScen5_modD[[i]] <- coxph (Surv (adj_surv_t, vital_st) ~ group + entry_t2,
                                 data = datasets5[[i]]$S_data, x = TRUE, model = TRUE,
                                 timefix = FALSE)
    
    # Longitudinal model on natural and adjusted time
    lmeScen5_time[[i]] <- lme(Y ~ time,
                              random = ~ time | ID, data = datasets5[[i]]$L_data,
                              control = lmeControl(opt = "optim"))
    lmeScen5_adjtime[[i]] <- lme(Y ~ adj_time,
                                 random = ~ adj_time | ID, data = datasets5[[i]]$L_data,
                                 control = lmeControl(opt = "optim"))
}

### Save Cox and LME models to file
saveRDS(coxScen5_modA, "Results/Sim_sub-models/coxScen5_modA.RDS")
saveRDS(coxScen5_modB, "Results/Sim_sub-models/coxScen5_modB.RDS")
saveRDS(coxScen5_modC, "Results/Sim_sub-models/coxScen5_modC.RDS")
saveRDS(coxScen5_modD, "Results/Sim_sub-models/coxScen5_modD.RDS")
saveRDS(lmeScen5_time, "Results/Sim_sub-models/lmeScen5_time.RDS")
saveRDS(lmeScen5_adjtime, "Results/Sim_sub-models/lmeScen5_adjtime.RDS")


# Note seed will be specified for ach run of future_lapply which manages the
# L'Ecuyer-CMRGbetween algorithm internally whne provided with a seed 

# Build joint models for scenario 5 model A all datasets
JM_Scen5_modA <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen5_time[[x]], coxScen5_modA[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 5623355)
# Save this batch
saveRDS(JM_Scen5_modA, "Results/Sim_Jms/JM_Scen5_modA.RDS")
# Clean-up to free up RAM
rm(JM_Scen5_modA)
gc()

###

# Build joint models for scenario 5 model B all datasets
JM_Scen5_modB <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen5_time[[x]], coxScen5_modB[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 83956244)
# Save this batch
saveRDS(JM_Scen5_modB, "Results/Sim_Jms/JM_Scen5_modB.RDS")
# Clean-up to free up RAM
rm(JM_Scen5_modB)
gc()

###

# Build joint models for scenario 5 model C all datasets
JM_Scen5_modC <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen5_adjtime[[x]], coxScen5_modC[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 35535)
# Save this batch
saveRDS(JM_Scen5_modC, "Results/Sim_Jms/JM_Scen5_modC.RDS")
# Clean-up to free up RAM
rm(JM_Scen5_modC)
gc()

###

# Build joint models for scenario 5 model D all datasets
JM_Scen5_modD <- future_lapply(1:M, function(x)
    jointModelBayes(lmeScen5_adjtime[[x]], coxScen5_modD[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 982956)
# Save this batch
saveRDS(JM_Scen5_modD, "Results/Sim_Jms/JM_Scen5_modD.RDS")
# Clean-up to free up RAM
rm(JM_Scen5_modD)
gc()



##### Build LME models with 2nd degree polynomials for non-linear JMs
poly_lme_time <- future_lapply(1:M, function(x)
    lme(Y ~ poly(time, 2),
            random=~ time | ID,
            data = datasets5[[x]]$L_data,
            control = lmeControl(opt = "optim")),
    future.seed = 34534)

poly_lme_adjtime <- future_lapply(1:M, function(x)
   lme(Y ~ poly(adj_time, 2),
       random = ~ adj_time | ID,
       data = datasets5[[x]]$L_data,
       control = lmeControl(opt = "optim")),
   future.seed = 7294575)

saveRDS(poly_lme_time, "Results/Sim_sub-models/poly_lme_time.RDS")
saveRDS(poly_lme_adjtime, "Results/Sim_sub-models/poly_lme_adjtime.RDS")


# Build joint models for scenario 5 model E all datasets
JM_Scen5_modE <- future_lapply(1:M, function(x)
    jointModelBayes(poly_lme_time[[x]], coxScen5_modB[[x]], timeVar = "time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 1920485)
# Save this batch
saveRDS(JM_Scen5_modE, "Results/Sim_Jms/JM_Scen5_modE.RDS")
# Clean-up to free up RAM
rm(JM_Scen5_modE)
gc()

# Build joint models for scenario 5 model F all datasets
JM_Scen5_modF <- future_lapply(1:M, function(x)
    jointModelBayes(poly_lme_adjtime[[x]], coxScen5_modD[[x]], timeVar = "adj_time",
                    control = list(n.iter = 10000, keepRE = FALSE)),
    future.seed = 1920485)
# Save this batch
saveRDS(JM_Scen5_modF, "Results/Sim_Jms/JM_Scen5_modF.RDS")
# Clean-up to free up RAM
rm(JM_Scen5_modF)
gc()



# Close down open multisession workers
plan(sequential)
