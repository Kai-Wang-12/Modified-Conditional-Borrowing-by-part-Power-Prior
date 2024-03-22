#########################################################################################################################################################################################################
#-----------------      This file includes R functions for leveraging external Gaussian endpoint with unknown variance from a single external control arm   -------------------#
#########################################################################################################################################################################################################
#    
# Following functions are included: 
#
# 1) Norm_convention() for leveraging external Gaussian endpoint with unknown variance by conventional Bayesian analysis;  
#
# 2) Norm_MPP() for leveraging external Gaussian endpoint with unknown variance by 1) Modified Power Prior with unknown (MPP1) and known variance (MPP2), 2) Modified Borrowing-by-part Power Prior (MBPP);      
#
# 3) Norm_MCBPP() for leveraging external Gaussian endpoint with unknown variance by the proposed Modified Conditional Borrowing-by-part Power Prior (MCBPP) and robust MCBPP (rMCBPP);   
#
# 4) Norm_CP_gamma() for leveraging external Gaussian endpoint with unknown variance by Commensurate Prior (CP) with conjugate gamma initial prior for the heterogeneity parameter;
#
# 5) Get_normal_rst1() for generating numerical results;
#
# 6) Get_normal_rst2() for generating simulation results;
#
# 7) Norm_numerical_function() for conducting the numerical study of leveraging external Gaussian endpoint with unknown variance from a single external control arm;   
#
# 8) Creat_dataset_normal() for generating simulated datasets for the simulation study of leveraging external Gaussian endpoint with unknown variance from a single external control arm;
#
# 9) Norm_simulation_function() for conducting the simulation study of leveraging external Gaussian endpoint with unknown variance from a single external control arm;   
#
# 10) Norm_application_function() for conducting the application.
#
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################r################################################################################################################################



############################################################################################
# Full Bayesian dynamic borrowing methods for leveraging external Gaussian endpoint with unknown variance from a single external control arm
############################################################################################

######### Inputs:

#@ N_CT: sample size of the experiment arm.
#@ N_CC: sample size of the concurrent controls.
#@ N_EC: sample size of the single external controls.
#
#@ Ybar_CT: sample mean of the experiment arm.
#@ Ybar_CT: sample mean of the concurrent controls.
#@ Ybar_EC: sample mean of the single external controls.

#@ Yvar_CT: sample variance of the experiment arm.
#@ Yvar_CC: sample variance of the concurrent controls.
#@ Yvar_EC: sample variance of the single external controls.
 
#@ type: type of full Bayesian dynamic borrowing methods, including 
#        a) "MPP1", "MPP2" for Norm_MPP(),
#        b) "MBPP" for Norm_MBPP(), and 
#        c) "MCBPP", "rMCBPP" for Norm_MCBPP().
#@ seed: provide value for the argument of "seed" in function stan() which generates initial values for MCMC.
#@ pars: provide value for the argument of "pars" in function stan() which specifies parameters of interest to be saved, where 
#        a) "mu_CT" denotes the population mean of the experiment arm, and 
#        b) "mu_CC", "sigma_CC" denotes the population mean and variance of the concurrent controls.
#@ chains: number of Markov chains (default = 1).
#@ alpha1, beta1: specify the conjugate gamma prior for the heterogeneity parameter of CP.

######### Outputs:
# 
# Normal_fit: the posterior results of the population mean and variance for a) the experimental arm, and b) the concurrent controls with no borrowing or pooling
# Normal_MPP_fit: the posterior results of the population mean for the concurrent controls by the MPP1, MPP2  
# Normal_MBPP_fit: the posterior results of the population mean for the concurrent controls by the MBPP  
# Normal_MCBPP_fit: the posterior results of the population mean for the concurrent controls by the MCBPP and rMCBPP  
# Normal_CP_gamma_fit: the posterior results of the population mean for the concurrent controls by the CP with conjugate gamma initial prior

Norm_convention <- function(N_CT, Ybar_CT, Yvar_CT, N_CC, Ybar_CC, Yvar_CC, pars = c("mu_CC","mu_CT","sigma_CC"), chains = 1, seed){
  
  data <- list(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_CT = N_CT, Ybar_CT = Ybar_CT, Yvar_CT = Yvar_CT)
  
  # 5000 iterations with a burn-in period of 500
  Normal_fit <- stan(
    file = "norm.stan", data = data, pars = pars, seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0
  )
  
  return(Normal_fit)
}
Norm_MPP        <- function(N_EC, N_CC, Ybar_EC, Ybar_CC, Yvar_CC, Yvar_EC, type, pars = c("mu_CC"), chains = 1, seed){
  
  data <- list(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC)
  
  if (type == "MPP1"){
    # 5000 iterations with a burn-in period of 500
    # delta: discounting parameter
    Normal_MPP_fit <- stan(file = "norm_MPP1.stan", data = data, pars = c(pars,"delta"),
                           seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  
  if (type == "MPP2"){
    # 5000 iterations with a burn-in period of 500
    # delta: discounting parameter
    Normal_MPP_fit <- stan(file = "norm_MPP2.stan", data = data, pars = c(pars,"delta"),
                           seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  
  return(Normal_MPP_fit)
  
}
Norm_MBPP        <- function(N_EC, N_CC, Ybar_EC, Ybar_CC, Yvar_CC, Yvar_EC, type, pars = c("mu_CC"), chains = 1, seed){
  
  data <- list(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC)
  if (type == "MBPP"){
    # 5000 iterations with a burn-in period of 500
    # delta_mu: mean-discounting parameter
    # delta_sigma: variance-discounting parameter
    Normal_MBPP_fit <- stan(file = "norm_MBPP.stan", data = data, pars = c(pars,"delta_mu","delta_sigma"),
                           seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  return(Normal_MBPP_fit)
  
}
Norm_MCBPP      <- function(N_EC, N_CC, Ybar_EC, Ybar_CC, Yvar_CC, Yvar_EC, type, pars = c("mu_CC"), chains = 1, seed){
  
  data <- list(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC)
  
  if (type == "MCBPP"){
    # 5000 iterations with a burn-in period of 500
    # delta_mu: mean-discounting parameter
    # delta_sigma: variance-discounting parameter
    Normal_MCBPP_fit <- stan(file = "norm_MCBPP.stan", data = data, pars = c(pars,"delta_mu","delta_sigma"),
                             seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  
  if (type == "rMCBPP"){
    # 5000 iterations with a burn-in period of 500
    # delta_mu: mean-discounting parameter
    Normal_MCBPP_fit <- stan(file = "norm_rMCBPP.stan", data = data, pars = c(pars,"delta_mu"),
                             seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }

  return(Normal_MCBPP_fit)
  
}
Norm_CP_gamma   <- function(N_EC, N_CC, Ybar_EC, Ybar_CC, Yvar_CC, Yvar_EC, type, pars = c("mu_CC"), alpha1 = 1, beta1 = 0.01, chains = 1, seed){
  
  data <- list(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, alpha1 = alpha1, beta1 = beta1)
  
  # 5000 iterations with a burn-in period of 500
  # tau2_mu: heterogeneity parameter
  Normal_CP_gamma_fit <- stan(file = "norm_CP_gamma.stan", data = data, pars = c(pars,"tau2_mu"),
                              seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  
  return(Normal_CP_gamma_fit)
  
}


############################################################################################
# Generating the results of numerical study 
############################################################################################

######### Inputs:

#@ post_C: posterior results of concurrent controls
#@ par_C: parameter name of the population mean for the concurrent controls
#@ par_mean: name of the mean-discounting parameter
#@ par_variance: name of variance-discounting parameter

######### Outputs:

# post_mean: the posterior mean of the mean-discounting parameter
# post_variance: the posterior mean of the variance-discounting parameter
# mean: the posterior mean of the population mean for the concurrent controls
# se: the posterior SE of the population mean for the concurrent controls

Get_normal_rst1 <- function(post_C, par_C = "mu_CC", par_mean = NA, par_variance = NA){
  
  data.frame(mean = mean(post_C[, par_C]), se = sd(post_C[, par_C]),
             post_mean = ifelse(is.na(par_mean), NA, mean(post_C[,par_mean])), 
             post_variance = ifelse(is.na(par_variance), NA, mean(post_C[,par_variance]))) %>%
  return()
  
}


############################################################################################
# Computing the results of the simulation study 
############################################################################################

######### Inputs:

#@ N_CC: sample size of the concurrent controls
#@ N_EC: sample size of the external controls
#@ post_T: posterior results of a) the population mean for the experiment arm or, b) the population mean and variance for the concurrent controls without borrowing
#@ par_T: parameter name of the population mean for the experiment arm in post_T
#@ par_C_scale: parameter name of the population variance for the concurrent controls without borrowing in post_T
#@ post_C: posterior results of the population mean for the concurrent controls with borrowing
#@ par_C: parameter name of the population mean for the concurrent controls in post_C
#@ norm_effect: the true value of average treatment effect

######### Outputs:

# ESS: the prior ESS
# Proportion: the proportion of prior ESS 
# Delta_mean: the posterior mean of the average treatment effect
# MSE: the posterior mse of the average treatment effect
# post_probability: the posterior probability of superiority
# CI_length: the length of 95% credible interval

Get_normal_rst2 <- function(N_CC, N_EC, post_T, par_T = "mu_CT", par_C_scale = "sigma_CC", post_C, par_C = "mu_CC", norm_effect = NA){
  
  set.seed(seed = 2023) # for the reproducibility of approximation
  scale <- sqrt(mean(post_T[,par_C_scale]))
  normmixd <- mixfit(post_C[, par_C], type = "norm", Nc = 2)
  ESS <- ess(normmixd, "elir", sigma = scale) - N_CC
  
  Delta  <- post_T[, par_T] - post_C[, par_C]
  
  post_probability  <- mean(Delta > 0)
  upper <- posterior_interval(matrix(Delta), prob = 0.95)[,2]
  lower <- posterior_interval(matrix(Delta), prob = 0.95)[,1]
  CI_length = upper - lower
  
  return(data.frame(ESS = ESS, Proportion = ESS/N_EC, Delta_mean = mean(Delta), 
                    MSE = ifelse(is.na(norm_effect),NA,mean((Delta-norm_effect)^2)),
                    post_probability = post_probability, CI_length = CI_length))
  
}


############################################################################################
# Conducting the numerical study of leveraging external Gaussian endpoint with unknown variance from a single external control arm
############################################################################################

######### Inputs:

#@ N: vector of the sample size of the experiment arm, concurrent controls and external controls
#@ Ybar: vector of the sample mean of the experiment arm, concurrent controls and external controls
#@ Yvar: vector of the sample variance of the experiment arm, concurrent controls and external controls
#@ seed: provide value for the argument of "seed" in function stan() which generates initial values for MCMC

######### Outputs:

# rst_n: the results of numerical study

Norm_numerical_function <- function(N, Ybar, Yvar, seed){
  
  N_CT = N[1];    N_CC = N[2];    N_EC = N[3]; 
  Ybar_CT = Ybar[1]; Ybar_CC = Ybar[2]; Ybar_EC = Ybar[3];
  Yvar_CT = Yvar[1]; Yvar_CC = Yvar[2]; Yvar_EC = Yvar[3];
  
  # 100 Markov chains 
  
  fit0_n <- Norm_convention(N_CT = N_CT, Ybar_CT = Ybar_CT, Yvar_CT = Yvar_CT, N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, chains = 100, seed = seed) %>% as.data.frame()
  fit1_n <- Norm_MPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MPP1", chains = 100, seed = seed) %>% as.data.frame()
  fit2_n <- Norm_MPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MPP2", chains = 100, seed = seed) %>% as.data.frame()
  fit3_n <- Norm_MBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MBPP", chains = 100, seed = seed) %>% as.data.frame()
  fit4_n <- Norm_MCBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MCBPP", chains = 100, seed = seed) %>% as.data.frame()
  fit5_n <- Norm_MCBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "rMCBPP", chains = 100, seed = seed) %>% as.data.frame()
  
  # 200 Markov chains 
  
  fit6_n <- Norm_CP_gamma(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, chains = 200, seed = seed) %>%as.data.frame()
  
  #
  rst0_n <- Get_normal_rst1(post_C = fit0_n) %>% data.frame(prior = "No-borrowing")
  rst1_n <- Get_normal_rst1(post_C = fit1_n, par_mean = "delta") %>% data.frame(prior = "MPP1")
  rst2_n <- Get_normal_rst1(post_C = fit2_n, par_mean = "delta") %>% data.frame(prior = "MPP2")
  rst3_n <- Get_normal_rst1(post_C = fit3_n, par_mean = "delta_mu", par_variance = "delta_sigma") %>% data.frame(prior = "MBPP")
  rst4_n <- Get_normal_rst1(post_C = fit4_n, par_mean = "delta_mu", par_variance = "delta_sigma") %>% data.frame(prior = "MCBPP")
  rst5_n <- Get_normal_rst1(post_C = fit5_n, par_mean = "delta_mu") %>% data.frame(prior = "rMCBPP")
  rst6_n <- Get_normal_rst1(post_C = fit6_n, par_mean = "tau2_mu") %>% data.frame(prior = "CP")
  
  rst_n <- rbind(rst0_n,rst1_n,rst2_n,rst3_n,rst4_n,rst5_n,rst6_n)
  return(rst_n)
}


############################################################################################
# Generating simulated datasets for the simulation study of leveraging external Gaussian endpoint with unknown variance from a single external control arm
############################################################################################

######### Inputs:

#@ N_CT: sample size of the experimental arm
#@ N_CC: sample size of the concurrent controls
#@ N_EC: sample size of the single external controls
#@ mu_CC: the population mean of the concurrent controls
#@ sd_C: equal population variance of the concurrent trial
#@ norm_effect: true average treatment effect
#@ LOC_hetero: vector of location heterogeneity
#@ VAR_hetero: the variance heterogeneity
#@ N_dataset: the number of simulated datasets
#@ seed: the seed for generating the simulated datasets

######### Outputs:

# data_norm: summary statistics of simulated datasets, including sample size, sample mean and variance of each arm and the pooled control arm.

Creat_dataset_SingleEC <- function(N_CC, N_CT, N_EC, mu_CC, LOC_hetero, norm_effect, sd_C, VAR_hetero, N_dataset, seed){
  
  data_norm <- list()
  
  ### Generate external controls with constraints
  
  set.seed(seed)
  
  Y_EC_div <- matrix(NA, nrow = length(LOC_hetero), ncol = N_EC)
  Ybar_EC_div <- rep(NA,length(LOC_hetero))
  Yvar_EC_div <- rep(NA,length(LOC_hetero))
  for (n in 1:length(LOC_hetero)){
    Y_EC <- lapply(as.list(1:10000),function(x){rnorm(N_EC, mu_CC + LOC_hetero[n]*sd_C, sd_C*VAR_hetero)})
    rst <- lapply(Y_EC, function(x){data.frame(m = abs(mean(x) - mu_CC - LOC_hetero[n]*sd_C), d = abs(sd(x) - sd_C*VAR_hetero))}) %>%
      rbindlist() %>%
      rowid_to_column() %>%
      mutate(sum = m+d) %>%
      filter(m < 0.05 & d < 0.1) %>%   # constraints
      arrange(sum) 
    Y_EC_div[n,] <- Y_EC[[as.numeric(rst[1,1])]]
    Ybar_EC_div[n] = mean(Y_EC_div[n,]); Yvar_EC_div[n] = var(Y_EC_div[n,])
    
  }
  
  ### Generate the concurrent trial
  
  set.seed(seed)
  
  for(i in 1:N_dataset){
    
    Y_CT1 <- rnorm(N_CT, mu_CC + norm_effect, sd_C)
    Y_CT0 <- rnorm(N_CT, mu_CC, sd_C)
    Y_CC <- rnorm(N_CC, mu_CC, sd_C)
    
    Ybar_CT1 = mean(Y_CT1); Yvar_CT1 = var(Y_CT1)
    Ybar_CT0 = mean(Y_CT0); Yvar_CT0 = var(Y_CT0)
    Ybar_CC = mean(Y_CC); Yvar_CC = var(Y_CC)
    
    data <- data.frame(N_CT = rep(N_CT,length(LOC_hetero)), N_CC = rep(N_CC,length(LOC_hetero)), N_EC = rep(N_EC,length(LOC_hetero)), 
                       Ybar_CT1 = rep(Ybar_CT1,length(LOC_hetero)), Ybar_CT0 = rep(Ybar_CT0,length(LOC_hetero)), 
                       Ybar_CC = rep(Ybar_CC,length(LOC_hetero)), Ybar_EC = rep(NA,length(LOC_hetero)), 
                       Yvar_CT1 = rep(Yvar_CT1,length(LOC_hetero)), Yvar_CT0 = rep(Yvar_CT0,length(LOC_hetero)),
                       Yvar_CC = rep(Yvar_CC,length(LOC_hetero)), Yvar_EC = rep(NA,length(LOC_hetero)), 
                       N_CC_all = rep(N_EC+N_CC,length(LOC_hetero)), Ybar_CC_all = rep(NA,length(LOC_hetero)), 
                       Yvar_CC_all = rep(NA,length(LOC_hetero)), heterogeneity = LOC_hetero)
    
    for (n in 1:length(LOC_hetero)){
      Ybar_CC_all = mean(c(Y_EC_div[n,],Y_CC)); Yvar_CC_all = var(c(Y_EC_div[n,],Y_CC))
      data[n,"Ybar_EC"] = Ybar_EC_div[n]; data[n,"Yvar_EC"] = Yvar_EC_div[n]
      data[n,"Ybar_CC_all"] = Ybar_CC_all; data[n,"Yvar_CC_all"] = Yvar_CC_all
    }
    
    data_norm[[i]] <- data
    
  }
  return(data_norm)
}


############################################################################################
# Conducting the simulation study of leveraging external Gaussian endpoint with unknown variance from a single external control arm
############################################################################################

######### Inputs:
#@ data: summary statistics of each simulated dataset
#@ Pooled: logical value indicating whether the "Pooled" should be conducted
#@ seed: provide value for the argument of "seed" in function stan() which generates initial values for MCMC

######### Outputs:
# rst_n: the results of each simulated dataset

Norm_simulation_function <- function(data, Pooled = T, seed){
  
  N_CT = data[1];    N_CC = data[2];    N_EC = data[3];   N_CC_all = data[12]
  Ybar_CT1 = data[4];Ybar_CT0 = data[5];Ybar_CC = data[6];Ybar_EC = data[7];Ybar_CC_all = data[13]
  Yvar_CT1 = data[8];Yvar_CT0 = data[9];Yvar_CC = data[10];Yvar_EC = data[11];Yvar_CC_all = data[14]
  
  ## analysis
  fit0_n_H0 <- Norm_convention(N_CT = N_CT, Ybar_CT = Ybar_CT0, Yvar_CT = Yvar_CT0, N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, seed = seed) %>% as.data.frame()
  fit0_n_H1 <- Norm_convention(N_CT = N_CT, Ybar_CT = Ybar_CT1, Yvar_CT = Yvar_CT1, N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, seed = seed) %>% as.data.frame()
  fit1_n <- Norm_MPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MPP1", seed = seed) %>% as.data.frame()
  fit2_n <- Norm_MPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MPP2", seed = seed) %>% as.data.frame()
  fit3_n <- Norm_MBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MBPP", seed = seed) %>% as.data.frame()
  fit4_n <- Norm_MCBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MCBPP", seed = seed) %>% as.data.frame()
  fit5_n <- Norm_MCBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "rMCBPP", seed = seed) %>% as.data.frame()
  fit6_n <- Norm_CP_gamma(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, seed = seed) %>% as.data.frame()
  
  if (Pooled==T){
    fit7_n_H0 <- Norm_convention(N_CC = N_CC_all, Ybar_CC = Ybar_CC_all, Yvar_CC = Yvar_CC_all, N_CT = N_CT, Ybar_CT = Ybar_CT0, Yvar_CT = Yvar_CT0, seed = seed) %>% as.data.frame()
    fit7_n_H1 <- Norm_convention(N_CC = N_CC_all, Ybar_CC = Ybar_CC_all, Yvar_CC = Yvar_CC_all, N_CT = N_CT, Ybar_CT = Ybar_CT1, Yvar_CT = Yvar_CT1, seed = seed) %>% as.data.frame()
    
  }
  
  ## summary-H0
  rst0_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H0, post_C = fit0_n_H0) %>% data.frame(prior = "No-borrowing", test = "H0")
  rst1_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H0, post_C = fit1_n) %>% data.frame(prior = "MPP1", test = "H0")
  rst2_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H0, post_C = fit2_n) %>% data.frame(prior = "MPP2", test = "H0")
  rst3_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H0, post_C = fit3_n) %>% data.frame(prior = "MBPP", test = "H0")
  rst4_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H0, post_C = fit4_n) %>% data.frame(prior = "MCBPP", test = "H0")
  rst5_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H0, post_C = fit5_n) %>% data.frame(prior = "rMCBPP", test = "H0")
  rst6_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H0, post_C = fit6_n) %>% data.frame(prior = "CP", test = "H0")
  
  
  ## summary-H1
  rst0_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H1, post_C = fit0_n_H1) %>% data.frame(prior = "No-borrowing", test = "H1")
  rst1_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H1, post_C = fit1_n) %>% data.frame(prior = "MPP1", test = "H1")
  rst2_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H1, post_C = fit2_n) %>% data.frame(prior = "MPP2", test = "H1")
  rst3_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H1, post_C = fit3_n) %>% data.frame(prior = "MBPP", test = "H1")
  rst4_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H1, post_C = fit4_n) %>% data.frame(prior = "MCBPP", test = "H1")
  rst5_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H1, post_C = fit5_n) %>% data.frame(prior = "rMCBPP", test = "H1")
  rst6_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H1, post_C = fit6_n) %>% data.frame(prior = "CP", test = "H1")
  
  if (Pooled==T){
    
    rst7_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H0, post_C = fit7_n_H0) %>% data.frame(prior = "Pooled", test = "H0")
    rst7_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n_H1, post_C = fit7_n_H1) %>% data.frame(prior = "Pooled", test = "H1")
    
    rst_n_H0 <- as.data.frame(rbind(rst0_n_H0,rst1_n_H0,rst2_n_H0,rst3_n_H0,rst4_n_H0,rst5_n_H0,rst6_n_H0,rst7_n_H0))
    rst_n_H1 <- as.data.frame(rbind(rst0_n_H1,rst1_n_H1,rst2_n_H1,rst3_n_H1,rst4_n_H1,rst5_n_H1,rst6_n_H1,rst7_n_H1))
    
  } else{
    rst_n_H0 <- as.data.frame(rbind(rst0_n_H0,rst1_n_H0,rst2_n_H0,rst3_n_H0,rst4_n_H0,rst5_n_H0,rst6_n_H0))
    rst_n_H1 <- as.data.frame(rbind(rst0_n_H1,rst1_n_H1,rst2_n_H1,rst3_n_H1,rst4_n_H1,rst5_n_H1,rst6_n_H1))
  }
  rst_n <- rbind(rst_n_H0, rst_n_H1)
  rownames(rst_n) <- NULL
  return(rst_n)
}


############################################################################################
# Generating simulated datasets for the application-----------------#
############################################################################################

######### Inputs:

#@ N_CT: sample size of the experimental arm
#@ N_CC: sample size of the concurrent controls
#@ N_EC: sample size of the single external controls
#@ mu_CC: the sample mean of the concurrent controls
#@ mu_EC: the sample mean of the external controls
#@ norm_effect: true average treatment effect (default = 0, H0)
#@ sd_CT: the sample standard deviation of the experimental arm
#@ sd_CC: the sample standard deviation of the concurrent controls
#@ sd_EC: the sample standard deviation of the external controls
#@ N_dataset: the number of simulated datasets
#@ seed: the seed for generating the simulated datasets

######### Outputs:

# data_norm: summary statistics of simulated datasets, including sample size, sample mean and variance.

Creat_dataset_normal_application <- function(N_CT, N_CC, N_EC, mu_EC, mu_CC, norm_effect = 0, sd_CT, sd_CC, sd_EC, N_dataset, seed){
  
  set.seed(seed)
  
  data <- data.frame(N_CT = rep(N_CT,N_dataset), N_CC = rep(N_CC,N_dataset), N_EC = rep(N_EC,N_dataset), 
                     Ybar_CT1 = rep(NA,N_dataset), Ybar_CT0 = rep(NA,N_dataset), Ybar_CC = rep(NA,N_dataset), Ybar_EC = rep(mu_EC,N_dataset), 
                     Yvar_CT1 = rep(NA,N_dataset), Yvar_CT0 = rep(NA,N_dataset), Yvar_CC = rep(NA,N_dataset), Yvar_EC = rep(sd_EC^2,N_dataset)) 
  
  for(i in 1:N_dataset){
    
    Y_CT1 <- rnorm(N_CT, mu_CC+norm_effect, sd_CT) 
    Y_CT0 <- rnorm(N_CT, mu_CC, sd_CT)
    Y_CC  <- rnorm(N_CC, mu_CC, sd_CC)
    
    Ybar_CT1 = mean(Y_CT1); Yvar_CT1 = var(Y_CT1)
    Ybar_CT0 = mean(Y_CT0); Yvar_CT0 = var(Y_CT0)
    Ybar_CC  = mean(Y_CC);  Yvar_CC  = var(Y_CC)
    data[i,"Ybar_CT1"] = Ybar_CT1; data[i,"Yvar_CT1"] = Yvar_CT1
    data[i,"Ybar_CT0"] = Ybar_CT0; data[i,"Yvar_CT0"] = Yvar_CT0
    data[i,"Ybar_CC"]  = Ybar_CC;  data[i,"Yvar_CC"]  = Yvar_CC
    
  }
  return(data)
}


############################################################################################
# Conducting the application
############################################################################################

######### Inputs:
#@ data: summary statistics of each simulated dataset
#@ norm_effect: the true value of average treatment effect
#@ chains: number of MCMC chains (default = 10)
#@ seed: provide value for the argument of "seed" in function stan() which generates initial values for MCMC

######### Outputs:
# rst_n: application results

Norm_application_function <- function(data, norm_effect, chains = 10, seed){
  
  N_CT = data[1];N_CC = data[2];N_EC = data[3]
  Ybar_CT = data[4];Ybar_CC = data[5];Ybar_EC = data[6]
  Yvar_CT = data[7];Yvar_CC = data[8];Yvar_EC = data[9]
  
  ## analysis
  fit0_n <- Norm_convention(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_CT = N_CT, Ybar_CT = Ybar_CT, Yvar_CT = Yvar_CT, chains = chains, seed = seed) %>% as.data.frame()
  fit1_n <- Norm_MPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MPP1", chains = chains, seed = seed) %>% as.data.frame()
  fit2_n <- Norm_MPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MPP2", chains = chains, seed = seed) %>% as.data.frame()
  fit3_n <- Norm_MBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MBPP", chains = chains, seed = seed) %>% as.data.frame()
  fit4_n <- Norm_MCBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "MCBPP", chains = chains, seed = seed) %>% as.data.frame()
  fit5_n <- Norm_MCBPP(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, type = "rMCBPP", chains = chains, seed = seed) %>% as.data.frame()
  fit6_n <- Norm_CP_gamma(N_EC = N_EC, N_CC = N_CC, Ybar_EC = Ybar_EC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, Yvar_EC = Yvar_EC, chains = chains, seed = seed) %>% as.data.frame()
  
  
  rst0_n <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n, post_C = fit0_n, norm_effect = norm_effect)
  rst1_n <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n, post_C = fit1_n, norm_effect = norm_effect)
  rst2_n <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n, post_C = fit2_n, norm_effect = norm_effect)
  rst3_n <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n, post_C = fit3_n, norm_effect = norm_effect)
  rst4_n <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n, post_C = fit4_n, norm_effect = norm_effect)
  rst5_n <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n, post_C = fit5_n, norm_effect = norm_effect)
  rst6_n <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC, post_T = fit0_n, post_C = fit6_n, norm_effect = norm_effect)
  
  rst_n <- as.data.frame(rbind(rst0_n,rst1_n,rst2_n,rst3_n,rst4_n,rst5_n,rst6_n))
  rst_n <- data.frame(prior = c("No-borrowing","MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP"), rst_n)
  rownames(rst_n) <- NULL
  return(rst_n)
}









