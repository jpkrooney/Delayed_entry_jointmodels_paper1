library(tidyverse)
library(mvtnorm)
library(JMbayes)
library(parallel)

source("RScripts/simJM_data.R")

# Set the number of dataset per scenario, M, to simulate
M <- 1000

# Set the number of cores to use for parallel computations
n_cpus <- 14

##### Set parameters for the simulations

# Number of patients in each subgroup
N <- c(500, 250)
# Weibull parameters
G0 <- -1.6 ; rho = 1.01
# Sub group increased hazards
group_Haz <- c(0, log(1.5))
# Left truncation parameters
ltrunc_meanlog <- 2.355
ltrunc_sdlog <- 0.768
# Association parameter between longitudinal and survival submodels
alpha <- -0.05
# Redisual variance
resid_var <- sqrt(2.4)
# Intercepts and slopes for longitudinal process
mu_lin = list( c (60, -1.0),
             c (60, -1.0))
mu_lin_gp = list( c (60, -1.0),
               c (60, -1.5))
mu_exp = list( c (60, 0.05),
               c (60, 0.05))
# Covariance matrix for longitudinal process
sig_vcov = matrix (c ( 20, -0.2,
                       -0.2, 0.25), 2)
vcov_exp = matrix (c ( 20, -0.002,
                       -0.002, 0.0002), 2)

# Follow up limit and max number of longitudinal values per individual
max_time <- 72 ; max_long_count = 10
# Paremeters to control longitudinal sampling frequency
interval <- 6 ; interval.sd = 2 ; interval.max = 48
#interval_nl <- 4; interval_nl_max = 60
# Hazard associated with delayed entry time
G_ltrunc_ninf <- 0    # non-informative
G_ltrunc_inf <- -0.05 # informative





##### Set RNG method and random seed
RNGkind(kind = "L'Ecuyer-CMRG")

set.seed(41673849)

##### Simulate each scenario

# Scenario 1
datasets1 <- mclapply(1: M, function(x)
    simJM_data(N = N, G0 = G0, rho = rho, group_Haz = group_Haz, alpha = alpha,
               resid_var = resid_var, mu_s = mu_lin, sig_vcov = sig_vcov,
               ltrunc = TRUE, G_ltrunc = G_ltrunc_ninf,
               ltrunc_meanlog = ltrunc_meanlog, ltrunc_sdlog = ltrunc_sdlog,
               max_time = max_time, interval = interval,
               interval.sd = interval.sd, interval.max = interval.max,
               cleanup = FALSE,
               long_mode = "linear", ltrunc_long_beta = 0),
    mc.cores = n_cpus, mc.set.seed = TRUE)

# Scenario 2
datasets2 <- mclapply(1: M, function(x)
    simJM_data(N = N, G0 = G0, rho = rho, group_Haz = group_Haz, alpha = alpha,
               resid_var = resid_var, mu_s = mu_lin, sig_vcov = sig_vcov,
               ltrunc = TRUE, G_ltrunc = G_ltrunc_inf,
               ltrunc_meanlog = ltrunc_meanlog, ltrunc_sdlog = ltrunc_sdlog,
               max_time = max_time, interval = interval,
               interval.sd = interval.sd, interval.max = interval.max,
               cleanup = FALSE,
               long_mode = "linear", ltrunc_long_beta = 0),
    mc.cores = n_cpus, mc.set.seed = TRUE)

# Scenario 3
datasets3 <- mclapply(1: M, function(x)
    simJM_data(N = N, G0 = G0, rho = rho, group_Haz = group_Haz, alpha = alpha,
               resid_var = resid_var, mu_s = mu_lin_gp, sig_vcov = sig_vcov,
               ltrunc = TRUE, G_ltrunc = G_ltrunc_inf,
               ltrunc_meanlog = ltrunc_meanlog, ltrunc_sdlog = ltrunc_sdlog,
               max_time = max_time, interval = interval,
               interval.sd = interval.sd, interval.max = interval.max,
               cleanup = FALSE,
               long_mode = "linear", ltrunc_long_beta = 0),
    mc.cores = n_cpus, mc.set.seed = TRUE)

# Scenario 4
datasets4 <- mclapply(1: M, function(x)
    simJM_data(N = N, G0 = G0, rho = rho, group_Haz = group_Haz, alpha = alpha,
               resid_var = resid_var, mu_s = mu_lin, sig_vcov = sig_vcov,
               ltrunc = TRUE, G_ltrunc = G_ltrunc_inf,
               ltrunc_meanlog = ltrunc_meanlog, ltrunc_sdlog = ltrunc_sdlog,
               max_time = max_time, interval = interval,
               interval.sd = interval.sd, interval.max = interval.max,
               cleanup = FALSE,
               long_mode = "linear", ltrunc_long_beta = 0.1),
    mc.cores = n_cpus, mc.set.seed = TRUE)

# Scenario 5
datasets5 <- mclapply(1: M, function(x)
    simJM_data(N = N, G0 = G0, rho = rho, group_Haz = group_Haz, alpha = alpha,
               resid_var = resid_var, mu_s = mu_exp, sig_vcov = vcov_exp,
               ltrunc = TRUE, G_ltrunc = G_ltrunc_inf,
               ltrunc_meanlog = ltrunc_meanlog, ltrunc_sdlog = ltrunc_sdlog,
               max_time = max_time, interval = interval,
               interval.sd = interval.sd, interval.max = interval.max,
               cleanup = FALSE,
               long_mode = "exp", ltrunc_long_beta = 0),
    mc.cores = n_cpus, mc.set.seed = TRUE)




##### Post simulation changes
# 1. Make a duplicate column of the dx_delay variable to faciliate modelling
for( i in 1:M){
    datasets1[[i]]$S_data$entry_t2 <- datasets1[[i]]$S_data$entry_t
    datasets1[[i]]$L_data$entry_t2 <- datasets1[[i]]$L_data$entry_t
    datasets2[[i]]$S_data$entry_t2 <- datasets2[[i]]$S_data$entry_t
    datasets2[[i]]$L_data$entry_t2 <- datasets2[[i]]$L_data$entry_t
    datasets3[[i]]$S_data$entry_t2 <- datasets3[[i]]$S_data$entry_t
    datasets3[[i]]$L_data$entry_t2 <- datasets3[[i]]$L_data$entry_t
    datasets4[[i]]$S_data$entry_t2 <- datasets4[[i]]$S_data$entry_t
    datasets4[[i]]$L_data$entry_t2 <- datasets4[[i]]$L_data$entry_t
    datasets5[[i]]$S_data$entry_t2 <- datasets5[[i]]$S_data$entry_t
    datasets5[[i]]$L_data$entry_t2 <- datasets5[[i]]$L_data$entry_t
}


# 2. Created adjusted time variable
for (i in 1:M){
    datasets1[[i]]$S_data$adj_surv_t <- datasets1[[i]]$S_data$surv_t - datasets1[[i]]$S_data$entry_t
    datasets1[[i]]$L_data$adj_time <- datasets1[[i]]$L_data$time - datasets1[[i]]$L_data$entry_t
    datasets2[[i]]$S_data$adj_surv_t <- datasets2[[i]]$S_data$surv_t - datasets2[[i]]$S_data$entry_t
    datasets2[[i]]$L_data$adj_time <- datasets2[[i]]$L_data$time - datasets2[[i]]$L_data$entry_t
    datasets3[[i]]$S_data$adj_surv_t <- datasets3[[i]]$S_data$surv_t - datasets3[[i]]$S_data$entry_t
    datasets3[[i]]$L_data$adj_time <- datasets3[[i]]$L_data$time - datasets3[[i]]$L_data$entry_t
    datasets4[[i]]$S_data$adj_surv_t <- datasets4[[i]]$S_data$surv_t - datasets4[[i]]$S_data$entry_t
    datasets4[[i]]$L_data$adj_time <- datasets4[[i]]$L_data$time - datasets4[[i]]$L_data$entry_t
    datasets5[[i]]$S_data$adj_surv_t <- datasets5[[i]]$S_data$surv_t - datasets5[[i]]$S_data$entry_t
    datasets5[[i]]$L_data$adj_time <- datasets5[[i]]$L_data$time - datasets5[[i]]$L_data$entry_t
    
}


# Save datasets
saveRDS(datasets1, "Data/sim_data1.RDS")
saveRDS(datasets2, "Data/sim_data2.RDS")
saveRDS(datasets3, "Data/sim_data3.RDS")
saveRDS(datasets4, "Data/sim_data4.RDS")
saveRDS(datasets5, "Data/sim_data5.RDS")




# Optionally visualise a longitudinal dataset
#ggplot(datasets5[[115]]$L_data, aes(x = time, y=Y, group=ID, col=group)) + geom_line()
## Optionally visualise a survival dataset
#plot (survfit (Surv (surv_t, vital_st) ~ group, data = datasets1[[5]]$S_data), col = 1:2,
#      ylab = "Survival probability", xlab = "Time since disease onset",
#      conf.int = FALSE)
#legend ("bottomleft", legend = c ("Group 1", "Group 2"), col = 1:2, lty = 1, bty = "n")




