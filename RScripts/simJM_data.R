#' Function to simulate joint longitudinal and time to event data. Default values are based on JM fit of longitudinal ALSFRS values and ALS survival data.
#' @param N vector defining number of sub-groups to simulate and number of individuals in each group.
#' @param G0 numeric.
#' @param group_Haz vector. Hazard associated with membership of a given group over baseline. group_Haz should be of same length as N. Negative number reduces event risk, positive increases event risk. Is equivalent to log(hazard ratio) for a given group. E.g. if simulating 2 groups where group 2 has HR = 2.0 with respect to group 1, log(2.0) = -0.4155154, therefore specify group_Haz as group_Haz = c(0, 0.6931472)
#' @param G_ltrunc numeric. Effect of delayed entry time on hazard.
#' @param max_time maximum time over which to simulate data.
#' @param interval mean time between each simulated longitudinal value.
#' @param interval.sd standard deviation of time between each simulated longtiduinal value.
#' @param interval.max last time at which to simulate longitudinal values.
#' @param max_long_count maximum number of longitudinal datapoints to simulate per individual. Default is 10.
#' @param ltrunc logical. If true data is simulated as left-truncated (delayed entry) data. This would typically be used to simulate observational data.
#' @param ltrunc_meanlog numeric. Left truncation times are simulated as lognormal variable. ltrunc_meanlog defines the mean on the log scale.
#' @param ltrunc_sdlog numeric. Left truncation times are simulated as lognormal variable. ltrunc_sdlog defines the standard deviation on the log scale.
#' @param ltrunc_long_beta numeric. Effect of left truncation time on the longitudinal process. Default = 0. If ltrunc = FALSE then ltrunc_long_beta has no effect.
#' @param group_long_factor numerical vector of same length as the number of subgroups defined in N. Default c(1, 1). To be deprecated - do not use. 
#' @param alpha numeric. Association of longitudinal and time-to-event processes.
#' @param rho numeric. Shape parameter of the Weibull distribution wihch defines the event process.
#' @param resid_var numeric.
#' @param mu_s list. For each sub-group in N, there should be a vector element of length 2 in mu_s. The first element of the vector defines the mean intercept of the longitudinal process with time = 0, the second defines the slope per unit time.
#' @param sig_vcov matrix. A 2x2 dimension matrix that gives variance covariance parameters for longitudinal intercept and slope simulation.
#' @param long_mode string. Default is "linear" which means the simulated longitudinal function has linear form. Other options are "exp" for an exponential decay.
#' @param censor_percent numeric. Percentage of simulated cohort to undergo censoring before death. Should be a number between 0 & 100. 10% be default.
#' @param cleanup boolean to indicate whether simulated data-sets should include the intermediate variables used to create the simulated data of not. Default is TRUE which means intermediate variables will not be returned.
#' @import dplyr mvtnorm magrittr
#' @export
# Main JM simulation function. Math and helper functions below.
# Based on code written by Ruben Van Eijk, University Medical Center Utrecht as part of the following paper: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5865572/
simJM_data <- function(N = c(200, 200),
                       G0 = -2.28, group_Haz = c(0, 0.6931472),
                       max_time = 36, max_long_count = 10,
                       interval = 3, interval.sd = 0.48, interval.max = 30,
                       ltrunc = FALSE, ltrunc_meanlog = NULL, ltrunc_sdlog = 0,
                       ltrunc_long_beta = 0, group_long_factor = c(1,1),
                       G_ltrunc = 0,
                       alpha = -0.12, rho = 1.44,
                       resid_var = sqrt (3.08),
                       mu_s = list( c (38, -1.06),
                                    c (38, -2.00)),
                       sig_vcov = matrix (c ( 26.10, 1.22,
                                              1.22, 1.28), 2),
                       long_mode = "linear",
                       censor_percent = 10,
                       cleanup = TRUE,
                       include_unobserved = FALSE){

    ##### Check arguments
    if( length(N) < 1 ){
        stop("Provide at least one N vector")
    }
    # add checks
    if (!long_mode %in% c("linear", "exp")) {
        stop("\'linear\' or \'exp\' are only valid options for variable \'long_mode\'.")
    }

    # To ensure consistent behaviour of group_long_factor by mode invert if mode = exp
    if(long_mode == "exp"){
        group_long_factor <- 1/ group_long_factor
    }

    ##### Generate survival data for groups
    all_groups <- list()
    # Make each group
    for(i in 1:length(N)){
        all_groups[[i]] <- make_survgroup(N, i=i, mu_s, sig_vcov, G0, group_Haz, G_ltrunc, alpha, rho,
                           ltrunc, ltrunc_meanlog, ltrunc_sdlog, ltrunc_long_beta, max_time, resid_var,
                           interval, interval.sd, interval.max, censor_percent, long_mode, group_long_factor)
    }
    surv_df <- do.call("rbind", all_groups)
    surv_df <- data.frame(ID = factor(1: nrow(surv_df)), surv_df)


    ##### Generate longitudinal data for surv_df
    # First generate time variable based on intervals + noise
    long_df <- data.frame(ID = rep (unique (surv_df$ID),
                                    each = length (seq (0, interval.max, by = interval))),
                          time = rep (seq (0, interval.max, by = interval), nrow(surv_df) ),
                          visit = rep (seq (0, interval.max, by = interval)/interval + 1, nrow(surv_df) ))
    long_df[long_df$time > 0, ]$time <- rnorm(nrow (long_df[long_df$time > 0, ]),
                                              mean = long_df[long_df$time > 0, ]$time,
                                              sd = interval.sd)
    # Create vital STATUS & survival TIME variable for survival analysis:
    surv_df$vital_st <- ifelse(surv_df$fail_t < surv_df$censor_t &
                                   surv_df$fail_t < max_time, 1, 0)
    surv_df$surv_t <- ifelse(surv_df$vital_st == 1, surv_df$fail_t, surv_df$censor_t)
    surv_df$surv_t <- ifelse(surv_df$surv_t > max_time, max_time, surv_df$surv_t)
    surv_df$group <- as.factor(surv_df$group)

    # Match baseline to longitudinal data
    long_df <- long_df %>% left_join( surv_df %>%
                                          dplyr::select( -c(fail_t, censor_t)),
                                      by = "ID")

    if(ltrunc == TRUE){
        # If observed == true for an individual make sure at least one longitudinal
        # time point exists within followup time
        # easiest way is to add one for every person
        temp <- long_df[ long_df$visit == 1 & long_df$observed == TRUE, ]
        # Check if any rows exist
        #if()
        temp$time <- runif( nrow(temp), min = temp$entry_t, max = temp$surv_t )
        temp$visit <- NA
        # join to main dataframe
        long_df <- rbind(long_df, temp)
        # sort by ID and time
        long_df <- long_df[ order(long_df$ID, long_df$time), ]
    }
    
    # Remove observed variable from longitudinal data as no longer needed and
    # may cause confusion with truncated longitudinal readings from individuals who are observed
    long_df$observed <- NULL
    
    
    ##### New generate longitudinal Y data at time t via appropriate equation
    if( long_mode == "linear"){
        long_df$Y <- sapply(1:nrow(long_df), function(x)
                    mi_lin( t1 = long_df[x, ]$time, B0 = long_df[x, ]$B0,
                        B1 = long_df[x, ]$B1, ltrunc_long_beta = long_df[x, ]$ltrunc_long_beta,
                        entry_t = long_df[x,]$entry_t)
                    ) + rnorm ( (nrow(long_df)), mean = 0, sd = long_df$resid_var)
    } else if (long_mode == "exp") {
        long_df$Y <- sapply(1:nrow(long_df), function(x)
            mi_exp( t1 = long_df[x, ]$time, B0 = long_df[x, ]$B0,
                    B1 = long_df[x, ]$B1, ltrunc_long_beta = long_df[x, ]$ltrunc_long_beta,
                    entry_t = long_df[x,]$entry_t)
        ) + rnorm ( (nrow(long_df)), mean = 0, sd = long_df$resid_var)
    }

    # Keep longitudinal values if survival time >= longitudinal time
    long_df <- long_df[ long_df$surv_t >= long_df$time, ]

    
    ###### If there is delayed entry some longitudinal variables will be unobserved
    # Will move these to a separate dataframe and remove from main dataframe
    if(ltrunc == TRUE) {
        trunc_long_df <- long_df[long_df$time < long_df$entry_t, ]
        long_df <- long_df[long_df$time >= long_df$entry_t, ]
        
        # Also remove unobserved individuals from survival dataframe into separate dataframe
        trunc_surv_df <- surv_df[surv_df$observed == FALSE, ]
        surv_df <- surv_df[surv_df$observed == TRUE, ]
    }
    
    # Reset visits counts to reflect additions or removals
    long_df <- long_df %>%
        group_by(ID) %>%
        mutate(visit = rank(time) )
    
    # Apply max_long_count to trim longitudinal data to max_long_count per individual
    long_df <- long_df[long_df$visit <= max_long_count, ]

    # If indvidual died at exact same time as last longitudinal follow up
    # remove small amount of time from the last longitudinal follow up
    # This avoids errors when building jointmodels
    long_df$time <- ifelse(long_df$time == long_df$surv_t,
                           long_df$surv_t - 0.0001, long_df$time)

    # Clean-up intermediate variables
    if( cleanup == TRUE){
        surv_df <- surv_df %>%
            dplyr::select( -c("B0", "B1", "resid_var", "G0", "group_Haz", "G_ltrunc",
                              "alpha", "rho", "s", "fail_t", "censor_t"))
        long_df <- long_df %>%
            dplyr::select( -c("B0", "B1", "resid_var", "G0", "group_Haz", "G_ltrunc",
                              "alpha", "rho", "s"))
    }

    # Package data for return to calling environment
    if((ltrunc == FALSE)){
        out <- list(S_data = data.frame(surv_df),
                    L_data = data.frame(long_df))
    } else {
    out <- list(S_data = data.frame(surv_df),
                L_data = data.frame(long_df),
                trunc_Sdata = data.frame(trunc_surv_df),
                trunc_Ldata = data.frame(trunc_long_df))
    }
    return (out)
}


# Function to make a population group
make_survgroup <- function(N, i, mu_s, sig_vcov, G0, group_Haz, G_ltrunc, alpha, rho,
                            ltrunc, ltrunc_meanlog, ltrunc_sdlog, ltrunc_long_beta, max_time, resid_var,
                           interval, interval.sd, interval.max, censor_percent, long_mode, group_long_factor){

    group <- list()
    # loop over each individual
    count = 0
    exclude_count = 0
    while(count < N[i] & exclude_count < N[i] ){
        ## Define longitudinal process vars for each individual:
        person <- as.data.frame (rmvnorm (1, mean = mu_s[[i]], sig_vcov))
        names(person) <- c ("B0", "B1")
        person$resid_var <- resid_var
        person$group <- i
        person$ltrunc_long_beta <- ltrunc_long_beta
        # Multiply B1 by group_long_factor depending in group
        person$group_long_factor <- group_long_factor[ person$group ]
        person$B1 <- person$B1 * person$group_long_factor

        # Add weibull parameters
        person$G0 <- G0; person$group_Haz <- group_Haz[i];  person$G_ltrunc <- G_ltrunc
        person$alpha <- alpha
        person$rho <- rho
        

        

        # Define a study entry time
        if(ltrunc == TRUE){
            person$entry_t <- rlnorm( 1, meanlog = ltrunc_meanlog, sdlog = ltrunc_sdlog)
        } else {
            person$entry_t <- 0
        }

        ### Calculate Time-to-event process, i.e. fail_t, by inverse transformation method
        person$s <- runif (1, min = 0, max = 1)
        # Search t1 between 0 and max_time months
        L <- try( uniroot(f = H, interval = c (0, max_time), S = person$s,
                          rho = person$rho, G0 = person$G0,
                          group_Haz = person$group_Haz, G_ltrunc = person$G_ltrunc,
                          ltrunc_t = person$entry_t,
                          alpha1 = person$alpha,
                          B0 = person$B0, B1 = person$B1,
                          ltrunc_long_beta = person$ltrunc_long_beta,
                          entry_t = person$entry_t,
                          long_mode )$root,
                  silent = TRUE)
        person$fail_t <- if (class (L) == "try-error"){
            max_time
        } else { L }

        ### Calculate a random censor time such that there is a 10% chance of being censored
        #person$censor_t <- runif (1) * rnorm (1, mean = 10*interval.max, sd = interval.sd)
        person$censor_t <- runif (1, 0, person$fail_t * (100/censor_percent) ) 
        
        
        ### Check if entry time occurs after censor_t or fail_t
        # if yes exclude this person
        if(person$entry_t > person$censor_t | person$entry_t > person$fail_t){
            person$observed <- FALSE
            exclude_count <- exclude_count + 1
            group[[count + exclude_count]] <- person
        } else {
            person$observed <- TRUE
            count = count + 1
            group[[count + exclude_count]] <- person
        }
    }
    group_df <- do.call("rbind", group)

    return(group_df)
}


### Math helper functions
# Linear longitudinal submodel
mi_lin <- function (t1, B0, B1, ltrunc_long_beta, entry_t){
    B0 + (B1 * t1) + (ltrunc_long_beta * entry_t)
}

# Exponential longitudinal submodel
mi_exp <- function(t1, B0, B1, ltrunc_long_beta, entry_t){
    (B0 * exp(-B1 * t1))  + (ltrunc_long_beta * entry_t)
}



# Inverse transformation method to calculate t1 for a given survival probability:
H <- function (t1, S, rho, G0, group_Haz, G_ltrunc, ltrunc_t, alpha1, B0, B1, long_mode,
               ltrunc_long_beta, entry_t){

    # Time-to-event submodel
    if(long_mode == "linear")
        h <- function (t1, rho, G0, group_Haz, G_ltrunc, ltrunc_t, alpha1,
                   B0, B1, ltrunc_long_beta, entry_t) {
            rho * exp (G0 + group_Haz + (G_ltrunc * ltrunc_t) +
                       (alpha1 * mi_lin(t1, B0, B1, ltrunc_long_beta, entry_t))) * (t1 ^ (rho - 1)
                       )
    } else if (long_mode == "exp") {
        h <- function (t1, rho, G0, group_Haz, G_ltrunc, ltrunc_t, alpha1,
                       B0, B1, ltrunc_long_beta, entry_t) {
            rho * exp (G0 + group_Haz + (G_ltrunc * ltrunc_t) +
                           (alpha1 * mi_exp(t1, B0, B1, ltrunc_long_beta, entry_t))) * (t1 ^ (rho - 1)
                           )
        }
    }
    # Calculate survival probability for time t1
    prob <- exp(-integrate(h, 0, upper=t1, rho=rho, G0=G0, group_Haz = group_Haz,
                            G_ltrunc = G_ltrunc, ltrunc_t=ltrunc_t, alpha1=alpha1,
                            B1=B1, B0=B0, ltrunc_long_beta, entry_t )$value) - S
    
    return(prob)
}

