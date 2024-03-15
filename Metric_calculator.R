# Some functions to quickly calculate RMSE, PBias, NSE, d, r, and r2 when multiple datasets in one 
# These functions were used in the manuscript: 
# "Model evaluation: the misuse of statistical techniques when evaluating observations versus predictions"
# Malcolm McPhee, Jonathan Richetti, Barry Croke, and Brad Walmsley
# 2024
# doi: xxx

# From hydroGOF (archived) https://cran.r-project.org/web/packages/hydroGOF/index.html
# Root Mean Square Error ----
rmse = function (sim, obs, na.rm = TRUE, ...) 
{
    # Root Mean Square Error - From hydroGOF (archived)
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    rmse <- sqrt(mean((obs - sim)^2, na.rm = na.rm))
    return(rmse)
}

# Percentage BIAS ----
pbias = function (sim, obs, na.rm = TRUE, ...) 
{
    # Percentage BIAS - From hydroGOF (archived)
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    n <- length(obs)
    denominator <- sum(obs)
    if (denominator != 0) {
        pbias <- 100 * (sum(sim - obs)/denominator)
    }
    else {
        pbias <- NA
        warning("'sum((obs)=0', it is not possible to compute 'pbias'")
    }
    return(round(pbias, 1))
}
valindex = function (sim, obs, ...) 
{
    if (length(obs) != length(sim)) {
        stop("Invalid argument: 'length(sim) != length(obs)' !! (", 
             length(sim), "!=", length(obs), ") !!")
    }
    else {
        index <- which(!is.na(sim) & !is.na(obs))
        if (length(index) == 0) 
            warning("'sim' and 'obs' are empty or they do not have any common pair of elements with data !!")
        return(index)
    }
}
# Nash-Sutcliffe efficiency ----
nse = function (sim, obs, na.rm = TRUE, FUN = NULL, epsilon = c(0, 
                                                                "Pushpalatha2012", "other"), epsilon.value = NA, ...) 
{
    #The Nash-Sutcliffe efficiency (NSE) is a normalized statistic that determines the relative magnitude of the residual variance ("noise") compared to the measured data variance ("information") (Nash and Sutcliffe, 1970).
    #Nash-Sutcliffe efficiencies range from -Inf to 1. Essentially, the closer to 1, the more accurate the model is.
    #-) NSE = 1, corresponds to a perfect match of modelled to the observed data.
    #-) NSE = 0, indicates that the model predictions are as accurate as the mean of the observed data,
    #-) -Inf < NSE < 0, indicates that the observed mean is better predictor than the model.
    #
    # From hydroGOF (archived)
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo", "xts"))) | is.na(match(class(obs), c("integer", 
                                                                              "numeric", "ts", "zoo", "xts")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo', 'xts')")
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.null(FUN)) {
        new <- preproc(sim = sim, obs = obs, FUN = FUN, epsilon = epsilon, 
                       epsilon.value = epsilon.value, ...)
        sim <- new[["sim"]]
        obs <- new[["obs"]]
    }
    denominator <- sum((obs - mean(obs))^2)
    if (denominator != 0) {
        NS <- 1 - (sum((obs - sim)^2)/denominator)
    }
    else {
        NS <- NA
        warning("'sum((obs - mean(obs))^2)=0' => it is not possible to compute 'NSE'")
    }
    return(NS)
}


# d-Willmott ----
d = function (sim, obs, na.rm = TRUE, ...) 
{
    #The Index of Agreement (d) developed by Willmott (1981) as a standardized measure of the degree of model prediction error and varies between 0 and 1.
    #A value of 1 indicates a perfect match, and 0 indicates no agreement at all (Willmott, 1981).
    #
    # From hydroGOF (archived)
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    Om <- mean(obs)
    denominator <- sum((abs(sim - Om) + abs(obs - Om))^2)
    if (denominator != 0) {
        d <- 1 - (sum((obs - sim)^2)/denominator)
    }
    else {
        d <- NA
        warning("'sum((abs(sim-Om)+abs(obs-Om))^2)=0', it is not possible to compute 'IoA'")
    }
    return(d)
}

# Pearsons r ----
rPearson = function (sim, obs, ...) 
{
    
    # Pearsons r
    # From hydroGOF (archived)
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    rPearson <- cor(sim, obs, method = "pearson", use = "pairwise.complete.obs")
    return(rPearson)
}

# 

#
# Modelling efficiency ----
mef = function (sim, obs, na.rm = TRUE, ...) 
{
    # MEF Modelling efficiency
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    mef = 1 - sum((obs - sim)^2)/sum((obs - mean(obs))^2)
    return(mef)
}

msep = function (sim, obs, na.rm = TRUE, ...) 
{
    # MSEP
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    msep = mean((obs - sim)^2)
    return(msep)
}

# Bias ----
msep_bias = function (sim, obs, na.rm = TRUE, ...) 
{
    # MSEP decomposed into error due to overall bias of prediction
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    msep_bias = (mean(sim) - mean(obs))^2
    return(msep_bias)
}

# Slope ----
msep_slope = function (sim, obs, na.rm = TRUE, ...) 
{
    # MSEP decomposed into error due to deviation of the regression slope from unity
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    
    ModFit = lm(obs ~ sim)
    SumFit = summary(ModFit)
    msep_slope = ((sum((sim - mean(sim))^2)) / length(obs)) * (1 - SumFit$coef[2])^2
    return(msep_slope)
}

#Deviance ----
msep_random = function (sim, obs, na.rm = TRUE, ...) 
{
    # MSEP decomposed into error due to overall bias of prediction
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    msep_random = (1 - cor(obs, sim)^2) * ((sum((obs - mean(obs))^2)) / length(obs))
    return(msep_random)
}


#Check sum of MSEP----
SumMSEP123 = function (sim, obs, na.rm = TRUE, ...) 
{
    # MSEP Check sum of MSEP
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    SumMSEP123 = msep_bias(obs, sim) + msep_slope(obs, sim) + msep_random(obs, sim)
    return(SumMSEP123)
}

# Bias proportion----
prop_bias = function (sim, obs, na.rm = TRUE, ...) 
{
    # Proportion of MSEP bias (%)
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    prop_bias = ( msep_bias(sim, obs) /  msep(sim, obs) ) *100
    return(prop_bias)
}

# Slope proportion----
prop_slope = function (sim, obs, na.rm = TRUE, ...) 
{
    # Proportion of MSEP Slope (%)
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    prop_slope = ( msep_slope(sim, obs) /  msep(sim, obs) ) *100
    return(prop_slope)
}

# Deviance proportion----
prop_deviance = function (sim, obs, na.rm = TRUE, ...) 
{
    # Proportion of MSEP Deviance (%)
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    prop_deviance = ( msep_random(sim, obs) /  msep(sim, obs) ) *100
    return(prop_deviance)
}


# T-test of the mean bias ----
P_mean_ttest = function (sim, obs, na.rm = TRUE, ...) 
{
    # returns the p-value of the Probability of paired t-test for the mean bias (P < 0.05).
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    P_slope_ttest = t.test(sim, obs, paired = TRUE)
    return(P_slope_ttest$p.value)
}

# beta coefficient  ----
beta = function (sim, obs, na.rm = TRUE, ...) 
{
    # returns the p-value of the Probability of paired t-test for the mean bias (P < 0.05).
    if (is.na(match(class(sim), c("integer", "numeric", "ts", 
                                  "zoo"))) | is.na(match(class(obs), c("integer", "numeric", 
                                                                       "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    if (length(obs) != length(sim)) 
        stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
    
    vi <- valindex(sim, obs)
    obs <- obs[vi]
    sim <- sim[vi]
    if (!is.na(match(class(sim), c("ts", "zoo")))) 
        sim <- as.numeric(sim)
    if (!is.na(match(class(obs), c("ts", "zoo")))) 
        obs <- as.numeric(obs)
    
    fitted_lm = summary(lm(obs ~ sim))
    beta = fitted_lm$coef[2]
    return(beta)
}

# wrapper functions ----

all_metrics_calc = function(df, obs, pred, group){
    # This function groups the dataframe in groups and then calculates the ALL metrics and return a data frame 
    library(tidyverse)

    group = enquo(group)
    obs = enquo(obs)
    pred = enquo(pred)
    
    # calculates the goodness of fit metrics for each group
    metrics = df %>%
        group_by(!!group) %>%
        reframe(n = length(!!pred),
                obs_mean = mean(!!obs),
                pred_mean = mean(!!pred),
                MeanBias = mean(!!obs - !!pred),
                MSEP = msep(!!obs, !!pred),
                RMSE = rmse(!!obs, !!pred),
                MSEP_bias = msep_bias(!!obs, !!pred),
                MSEP_slope = msep_slope(!!obs, !!pred),
                MSEP_random = msep_random(!!obs, !!pred),
                
                SumMSEP = SumMSEP123(!!obs, !!pred),
                Prop_BIAS = prop_bias(!!obs, !!pred),
                Prop_SLOPE = prop_slope(!!obs, !!pred),
                Prop_DEVIANCE = prop_deviance(!!obs, !!pred),
                
                MEF = mef(!!obs, !!pred),

                PBIAS = pbias(!!obs, !!pred),
                d = d(!!obs, !!pred),
                NSE = nse(!!obs, !!pred),
                r = rPearson(!!obs, !!pred),
                r2 = r^2,
                P_mean = P_mean_ttest(!!obs, !!pred),
                beta_coef = beta(!!obs, !!pred)
                )
    
    # a bit of fluffing around to get it a bit nicer :P
    transposed_metrics = t(metrics)
    transposed_metrics = transposed_metrics[-1,]
    
    # return the goodness of fit metrics
    return(transposed_metrics)
}

recomended_metrics_calc = function(df, obs, pred, group){
    # This function groups the dataframe in groups and then calculates the ALL metrics and return a data frame 
    library(tidyverse)
    
    group = enquo(group)
    obs = enquo(obs)
    pred = enquo(pred)
    
    # calculates the goodness of fit metrics for each group
    metrics = df %>%
        group_by(!!group) %>%
        reframe(n = length(!!pred),
                obs_mean = mean(!!obs),
                pred_mean = mean(!!pred),
                MeanBias = mean(!!obs - !!pred),
                RMSE = rmse(!!pred, !!obs),
                Prop_BIAS = prop_bias(!!pred, !!obs),
                Prop_SLOPE = prop_slope(!!pred, !!obs),
                Prop_DEVIANCE = prop_deviance(!!pred, !!obs),
                NSE_MEF = mef(!!pred, !!obs),
        )
    
    # a bit of fluffing around to get it a bit nicer :P
    transposed_metrics = t(metrics)
    transposed_metrics = transposed_metrics[-1,]
    
    # return the goodness of fit metrics
    return(transposed_metrics)
}

metrics_calc = function(df, obs, pred, group){
    # This function groups the dataframe in groups and then calculates the metrics and return a data frame 
    library(tidyverse)
    
    group = enquo(group)
    obs = enquo(obs)
    pred = enquo(pred)
    
    # calculates the goodness of fit metrics for each group
    metrics = df %>%
        group_by(!!group) %>%
        reframe(n = length(!!pred),
                RMSE = rmse(!!obs, !!pred),
                PBIAS = pbias(!!obs, !!pred),
                d = d(!!obs, !!pred),
                NSE = nse(!!obs, !!pred),
                r2 = rPearson(!!obs, !!pred)^2
        )
    
    # a bit of fluffing around to get it a bit nicer :P
    transposed_metrics = t(metrics)
    transposed_metrics = transposed_metrics[-1,]
    
    # return the goodness of fit metrics
    return(transposed_metrics)
}
#########################################################
## little example on how the function expects the data ##

#datasets = rep(1:4, each=5)
#obs = runif(n=20, min=1, max=20)
#pred = runif(n=20, min=1, max=20)
#df = data.frame(datasets,obs,pred)
#names(df) = c('Datasets','Observations','Simulations')

#all_metrics_calc(df, Observations, Simulations, Datasets)
