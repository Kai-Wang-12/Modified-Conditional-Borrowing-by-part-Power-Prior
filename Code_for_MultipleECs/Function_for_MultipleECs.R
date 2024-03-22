#########################################################################################################################################################################################################
#-------------     This file includes R functions for leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms   ---------------#
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
# 7) Norm_multiple_numerical_function() for conducting the numerical study of leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms;   
#
# 8) Creat_dataset_MultipleECs() for generating simulated datasets for the simulation study of leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms;
#
# 9) Norm_multiple_simulation_function() for conducting the simulation study of leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms.   
#
#########################################################################################################################################################################################################
#########################################################################################################################################################################################################
#########################################################################r################################################################################################################################



############################################################################################
# Full Bayesian dynamic borrowing methods for leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms
############################################################################################

######### Inputs:

#@ N_CT: sample size of the experimental arm.
#@ N_CC: sample size of the concurrent controls.
#@ N_EC1, N_EC2, N_EC3: sample size of the external control arm 1, 2 and 3.
#
#@ Ybar_CT: sample mean of the experimental arm.
#@ Ybar_CT: sample mean of the concurrent controls.
#@ Ybar_EC1, Ybar_EC2, Ybar_EC3: sample mean of the external control arm 1, 2 and 3. 

#@ Yvar_CT: sample variance of the experimental arm.
#@ Yvar_CC: sample variance of the concurrent controls.
#@ Yvar_EC1, Yvar_EC2, Yvar_EC3: sample variance of the external control arm 1, 2 and 3. 

#@ type: type of full Bayesian dynamic borrowing methods, including 
#        a) "MPP1_I", "MPP2_I" under Integrated Strategy and "MPP1_A", "MPP2_A" under Arm-based Strategy for Norm_MPP(),
#        b) "MBPP_I" under Integrated Strategy and "MBPP_A" under Arm-based Strategy for Norm_MBPP(), and 
#        c) "MCBPP_I" under Integrated Strategy and "MCBPP_A", "rMCBPP_A" under Arm-based Strategy for Norm_MCBPP().
#@ seed: provide value for the argument of "seed" in function stan() which generates initial values for MCMC
#@ pars: provide value for the argument of "pars" in function stan() which specifies parameters of interest to be saved, where 
#        a) "mu_CT" denotes the population mean of the experiment arm, and 
#        b) "mu_CC", "sigma_CC" denote the population mean and variance of the concurrent controls.
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
Norm_multi_MPP <- function(N_CC, Ybar_CC, Yvar_CC, N_EC1, Ybar_EC1, Yvar_EC1, N_EC2, Ybar_EC2, Yvar_EC2, N_EC3, Ybar_EC3, Yvar_EC3, 
                           pars = "mu_CC", chains = 1, type, seed){
  
  data <- list(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
               N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3)
  
  if (type == "MPP1_I"){
    # 5000 iterations with a burn-in period of 500
    # delta1, delta2, delta3: discounting parameter specific to EC1, EC2, EC3 
    Normal_MPP_fit <- stan(file = "norm_MPP1_I.stan", data = data, pars = c(pars,"delta1", "delta2", "delta3"),
                           seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  if (type == "MPP1_A"){
    # 5000 iterations with a burn-in period of 500
    # delta1, delta2, delta3: discounting parameter specific to EC1, EC2, EC3 
    Normal_MPP_fit <- stan(file = "norm_MPP1_A.stan", data = data, pars = c(pars,"delta1", "delta2", "delta3"),
                           seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  
  if (type == "MPP2_I"){
    # 5000 iterations with a burn-in period of 500
    # delta1, delta2, delta3: discounting parameter specific to EC1, EC2, EC3 
    Normal_MPP_fit <- stan(file = "norm_MPP2_I.stan", data = data, pars = c(pars,"delta1", "delta2", "delta3"),
                           seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  if (type == "MPP2_A"){
    # 5000 iterations with a burn-in period of 500
    # delta1, delta2, delta3: discounting parameter specific to EC1, EC2, EC3 
    Normal_MPP_fit <- stan(file = "norm_MPP2_A.stan", data = data, pars = c(pars,"delta1", "delta2", "delta3"),
                           seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0) 
  }
  
  return(Normal_MPP_fit)
  
}
Norm_multi_MBPP <- function(N_CC, Ybar_CC, Yvar_CC, N_EC1, Ybar_EC1, Yvar_EC1, N_EC2, Ybar_EC2, Yvar_EC2, N_EC3, Ybar_EC3, Yvar_EC3, 
                            pars = "mu_CC", chains = 1, type, seed){
  
  data <- list(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
               N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3)
  
  if (type == "MBPP_I"){
    # 5000 iterations with a burn-in period of 500
    # delta_mu1, delta_mu2, delta_mu3: the mean-discounting parameter specific to EC1, EC2, EC3 
    # delta_sigma1, delta_sigma2, delta_sigma3: the variance-discounting parameter specific to EC1, EC2, EC3 
    Normal_MBPP_fit <- stan(file = "norm_MBPP_I.stan", data = data, pars = c(pars,"delta_mu1", "delta_sigma1", "delta_mu2", "delta_sigma2", "delta_mu3", "delta_sigma3"),
                            seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  if (type == "MBPP_A"){
    # 5000 iterations with a burn-in period of 500
    # delta_mu1, delta_mu2, delta_mu3: the mean-discounting parameter specific to EC1, EC2, EC3 
    # delta_sigma1, delta_sigma2, delta_sigma3: the variance-discounting parameter specific to EC1, EC2, EC3 
    Normal_MBPP_fit <- stan(file = "norm_MBPP_A.stan", data = data, pars = c(pars,"delta_mu1", "delta_sigma1", "delta_mu2", "delta_sigma2", "delta_mu3", "delta_sigma3"),
                            seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0)
  }
  
  return(Normal_MBPP_fit)
  
}
Norm_multi_MCBPP<- function(N_CC, Ybar_CC, Yvar_CC, N_EC1, Ybar_EC1, Yvar_EC1, N_EC2, Ybar_EC2, Yvar_EC2, N_EC3, Ybar_EC3, Yvar_EC3, 
                            pars = "mu_CC", chains = 1, type, seed){
  
  data <- list(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
               N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3)
  
  if (type == "MCBPP_I"){
    # 5000 iterations with a burn-in period of 500
    # delta_mu1, delta_mu2, delta_mu3: the mean-discounting parameter specific to EC1, EC2, EC3 
    # delta_sigma1, delta_sigma2, delta_sigma3: the variance-discounting parameter specific to EC1, EC2, EC3  
    Normal_MCBPP_fit <- stan(file = "norm_MCBPP_I.stan", data = data, pars = c(pars,"delta_mu1", "delta_sigma1", "delta_mu2", "delta_sigma2", "delta_mu3", "delta_sigma3"),
                             seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0) 
  }
  if (type == "MCBPP_A"){
    # 5000 iterations with a burn-in period of 500
    # delta_mu1, delta_mu2, delta_mu3: the mean-discounting parameter specific to EC1, EC2, EC3  
    # delta_sigma1, delta_sigma2, delta_sigma3: the variance-discounting parameter specific to EC1, EC2, EC3 
    Normal_MCBPP_fit <- stan(file = "norm_MCBPP_A.stan", data = data, pars = c(pars,"delta_mu1", "delta_sigma1", "delta_mu2", "delta_sigma2", "delta_mu3", "delta_sigma3"),
                             seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0) 
  }
  
  if (type == "rMCBPP_A"){
    # 5000 iterations with a burn-in period of 500
    # delta_mu1, delta_mu2, delta_mu3: the mean-discounting parameter specific to EC1, EC2, EC3 
    Normal_MCBPP_fit <- stan(file = "norm_rMCBPP_A.stan", data = data, pars = c(pars,"delta_mu1", "delta_mu2", "delta_mu3"),
                             seed = seed, chains = chains, warmup = 500, iter = 5000, cores = 1, refresh = 0) 
  }
  
  return(Normal_MCBPP_fit)
  
}
Norm_multi_CP_gamma <- function(N_CC, Ybar_CC, Yvar_CC, N_EC1, Ybar_EC1, Yvar_EC1, N_EC2, Ybar_EC2, Yvar_EC2, N_EC3, Ybar_EC3, Yvar_EC3, 
                                pars = "mu_CC",  alpha1 = 1, beta1 = 0.01, chains = 1, seed){
  
  data <- list(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
               N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3,
               alpha1 = alpha1, beta1 = beta1)
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

# post_mean1, post_mean2, post_mean3: the posterior mean of mean-discounting parameter specific to EC1, EC2, EC3 
# post_variance1, post_variance2, post_variance3: the posterior mean of scale-discounting parameter specific to EC1, EC2, EC3 
# mean: the posterior mean of the population mean for concurrent controls
# se: the posterior SE of the population variance for concurrent controls

Get_normal_rst1 <- function(post_C, par_C = "mu_CC", par_mean = rep(NA,3), par_variance = rep(NA,3)){
  
  data.frame(mean = mean(post_C[, par_C]), se = sd(post_C[, par_C]),
             post_mean1 = ifelse(is.na(par_mean[1]), NA, mean(post_C[,par_mean[1]])), 
             post_variance1 = ifelse(is.na(par_variance[1]), NA, mean(post_C[,par_variance[1]])),
             post_mean2 = ifelse(is.na(par_mean[2]), NA, mean(post_C[,par_mean[2]])), 
             post_variance2 = ifelse(is.na(par_variance[2]), NA, mean(post_C[,par_variance[2]])),
             post_mean3 = ifelse(is.na(par_mean[3]), NA, mean(post_C[,par_mean[3]])), 
             post_variance3 = ifelse(is.na(par_variance[3]), NA, mean(post_C[,par_variance[3]]))) %>%
    return()
  
}


############################################################################################
# Generating the results of the simulation study 
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
  
  set.seed(seed = 2013) # for the reproducibility of the approximation
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
# Conducting the numerical study of leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms
############################################################################################

######### Inputs:

#@ N: vector of the sample size of the experiment arm, concurrent controls and external controls
#@ Ybar: vector of the sample mean of the experiment arm, concurrent controls and external controls
#@ Yvar: vector of the sample variance of the experiment arm, concurrent controls and external controls
#@ seed: provide value for the argument of "seed" in function stan() which generates initial values for MCMC

######### Outputs:

# rst_n: the results of numerical study

Norm_multiple_numerical_function <- function(N, Ybar, Yvar, seed){
  
  N_CT = N[1];    N_CC = N[2];    N_EC1 = N[3];    N_EC2 = N[4];    N_EC3 = N[5]; 
  Ybar_CT = Ybar[1]; Ybar_CC = Ybar[2]; Ybar_EC1 = Ybar[3];  Ybar_EC2 = Ybar[4];  Ybar_EC3 = Ybar[5];
  Yvar_CT = Yvar[1]; Yvar_CC = Yvar[2]; Yvar_EC1 = Yvar[3];  Yvar_EC2 = Yvar[4];  Yvar_EC3 = Yvar[5];
  
  # 100 Markov chains 
  
  fit0_n   <- Norm_convention(N_CT = N_CT, Ybar_CT = Ybar_CT, Yvar_CT = Yvar_CT, N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, chains = 100, seed = seed) %>% as.data.frame()
  fit1_n_I <- Norm_multi_MPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                             N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MPP1_I", chains = 100, seed = seed) %>% as.data.frame()
  fit1_n_A <- Norm_multi_MPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                             N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MPP1_A", chains = 100, seed = seed) %>% as.data.frame() 
  fit2_n_I <- Norm_multi_MPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                             N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MPP2_I", chains = 100, seed = seed) %>% as.data.frame() 
  fit2_n_A <- Norm_multi_MPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                             N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MPP2_A", chains = 100, seed = seed) %>% as.data.frame()
  fit3_n_I <- Norm_multi_MBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                              N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MBPP_I", chains = 100, seed = seed) %>% as.data.frame() 
  fit3_n_A <- Norm_multi_MBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                              N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MBPP_A", chains = 100, seed = seed) %>% as.data.frame() 
  fit4_n_I <- Norm_multi_MCBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                               N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MCBPP_I", chains = 100, seed = seed)  %>% as.data.frame()
  fit4_n_A <- Norm_multi_MCBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                               N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MCBPP_A", chains = 100, seed = seed) %>% as.data.frame()
  fit5_n_A <- Norm_multi_MCBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                               N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "rMCBPP_A", chains = 100, seed = seed) %>% as.data.frame()   
  
  # 200 Markov chains 
  
  fit6_n <- Norm_multi_CP_gamma(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                                N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, chains = 200, seed = seed)  %>% as.data.frame()
  
  
  rst0_n    <- Get_normal_rst1(post_C = fit0_n) %>% data.frame(prior = "No-borrowing")
  rst1_n_I  <- Get_normal_rst1(post_C = fit1_n_I,par_mean = c("delta1","delta2","delta3")) %>% data.frame(prior = "MPP1_I")
  rst1_n_A  <- Get_normal_rst1(post_C = fit1_n_A,par_mean = c("delta1","delta2","delta3")) %>% data.frame(prior = "MPP1_A")
  rst2_n_I  <- Get_normal_rst1(post_C = fit2_n_I,par_mean = c("delta1","delta2","delta3")) %>% data.frame(prior = "MPP2_I")
  rst2_n_A  <- Get_normal_rst1(post_C = fit2_n_A,par_mean = c("delta1","delta2","delta3")) %>% data.frame(prior = "MPP2_A")
  rst3_n_I  <- Get_normal_rst1(post_C = fit3_n_I,par_mean = c("delta_mu1","delta_mu2","delta_mu3"), 
                               par_variance = c("delta_sigma1","delta_sigma2","delta_sigma3")) %>% data.frame(prior = "MBPP_I")
  rst3_n_A  <- Get_normal_rst1(post_C = fit3_n_A,par_mean = c("delta_mu1","delta_mu2","delta_mu3"), 
                               par_variance = c("delta_sigma1","delta_sigma2","delta_sigma3")) %>% data.frame(prior = "MBPP_A")
  rst4_n_I  <- Get_normal_rst1(post_C = fit4_n_I,par_mean = c("delta_mu1","delta_mu2","delta_mu3"), 
                               par_variance = c("delta_sigma1","delta_sigma2","delta_sigma3")) %>% data.frame(prior = "MCBPP_I")
  rst4_n_A  <- Get_normal_rst1(post_C = fit4_n_A,par_mean = c("delta_mu1","delta_mu2","delta_mu3"), 
                               par_variance = c("delta_sigma1","delta_sigma2","delta_sigma3")) %>% data.frame(prior = "MCBPP_A")
  rst5_n_A  <- Get_normal_rst1(post_C = fit5_n_A,par_mean = c("delta_mu1","delta_mu2","delta_mu3")) %>% data.frame(prior = "rMCBPP_A")
  rst6_n    <- Get_normal_rst1(post_C = fit6_n,par_mean = "tau2_mu") %>% data.frame(prior = "CP")
  
  rst_n1 <- rbind(rst0_n,rst1_n_I,rst2_n_I,rst3_n_I,rst4_n_I,rst6_n) %>% mutate(method = "I")
  rst_n2 <- rbind(rst0_n,rst1_n_A,rst2_n_A,rst3_n_A,rst4_n_A,rst5_n_A) %>% mutate(method = "A")
  return(rbind(rst_n1,rst_n2))
}


############################################################################################
# Generating simulated datasets for the simulation study of leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms
############################################################################################

######### Inputs:

#@ N_CT: sample size of the experimental arm
#@ N_CC: sample size of the concurrent controls
#@ N_EC1, N_EC2, N_EC3: sample size of the external control arm 1, 2 and 3 
#@ mu_CC: the population mean of the concurrent controls
#@ mu_EC1,mu_EC2,mu_EC3: baseline value of the population mean for the external control arm 1, 2 and 3 
#@ sd_C: equal population variance of the concurrent trial
#@ norm_effect: the true value of average treatment effect
#@ sd_EC1, sd_EC2, sd_EC3: the population variance of the external control arm 1, 2 and 3 
#@ LOC_hetero: the vector of location heterogeneity
#@ N_dataset: the number of simulated datasets
#@ seed: the seed for generating the simulated datasets

######### Outputs:

# data_norm: summary statistics of simulated datasets, including sample size, sample mean and variance of each arm and the pooled control arm.

Creat_dataset_MultipleECs <- function(N_CC, N_CT, N_EC1, N_EC2, N_EC3, 
                                      mu_CC, mu_EC1, mu_EC2, mu_EC3, LOC_hetero, norm_effect, 
                                      sd_C, sd_EC1, sd_EC2, sd_EC3, N_dataset, seed){
  
  data_norm <- list()
  
  ### Generate external controls with constraints
  
  set.seed(seed)
  
  Y_EC_div1 <- matrix(NA, nrow = length(LOC_hetero), ncol = N_EC1)
  Ybar_EC_div1 <- rep(NA,length(LOC_hetero))
  Yvar_EC_div1 <- rep(NA,length(LOC_hetero))
  Y_EC_div2 <- matrix(NA, nrow = length(LOC_hetero), ncol = N_EC2)
  Ybar_EC_div2 <- rep(NA,length(LOC_hetero))
  Yvar_EC_div2 <- rep(NA,length(LOC_hetero))
  Y_EC_div3 <- matrix(NA, nrow = length(LOC_hetero), ncol = N_EC3)
  Ybar_EC_div3 <- rep(NA,length(LOC_hetero))
  Yvar_EC_div3 <- rep(NA,length(LOC_hetero))
  
  for (n in 1:length(LOC_hetero)){
    Y_EC1 <- lapply(as.list(1:10000),function(x){rnorm(N_EC1, mu_EC1 + LOC_hetero[n]*sd_C, sd_EC1)})
    rst1 <- lapply(Y_EC1, function(x){data.frame(m = abs(mean(x) - mu_EC1 - LOC_hetero[n]*sd_C), d = abs(sd(x) - sd_EC1))}) %>%
      rbindlist() %>%
      rowid_to_column() %>%
      mutate(sum = m+d) %>%
      filter(m < 0.05 & d < 0.05) %>%
      arrange(sum) 
    Y_EC_div1[n,] <- Y_EC1[[as.numeric(rst1[1,1])]]
    Ybar_EC_div1[n] = mean(Y_EC_div1[n,])
    Yvar_EC_div1[n] = var(Y_EC_div1[n,])
    
    Y_EC2 <- lapply(as.list(1:10000),function(x){rnorm(N_EC2, mu_EC2 + LOC_hetero[n]*sd_C, sd_EC2)})
    rst2 <- lapply(Y_EC2, function(x){data.frame(m = abs(mean(x) - mu_EC2 - LOC_hetero[n]*sd_C), d = abs(sd(x) - sd_EC2))}) %>%
      rbindlist() %>%
      rowid_to_column() %>%
      mutate(sum = m+d) %>%
      filter(m < 0.05 & d < 0.05) %>%
      arrange(sum) 
    Y_EC_div2[n,] <- Y_EC2[[as.numeric(rst2[1,1])]]
    Ybar_EC_div2[n] = mean(Y_EC_div2[n,])
    Yvar_EC_div2[n] = var(Y_EC_div2[n,])
    
    Y_EC3 <- lapply(as.list(1:10000),function(x){rnorm(N_EC3, mu_EC3 + LOC_hetero[n]*sd_C, sd_EC3)})
    rst3 <- lapply(Y_EC3, function(x){data.frame(m = abs(mean(x) - mu_EC3 - LOC_hetero[n]*sd_C), d = abs(sd(x) - sd_EC3))}) %>%
      rbindlist() %>%
      rowid_to_column() %>%
      mutate(sum = m+d) %>%
      filter(m < 0.05 & d < 0.05) %>%
      arrange(sum) 
    Y_EC_div3[n,] <- Y_EC3[[as.numeric(rst3[1,1])]]
    Ybar_EC_div3[n] = mean(Y_EC_div3[n,])
    Yvar_EC_div3[n] = var(Y_EC_div3[n,])
  }
  
  
  ### Generate concurrent trial
  
  for(i in 1:N_dataset){
    
    Y_CT1 <- rnorm(N_CT, mu_CC + norm_effect, sd_C)
    Ybar_CT1 = mean(Y_CT1); Yvar_CT1 = var(Y_CT1)
    
    Y_CT0 <- rnorm(N_CT, mu_CC, sd_C)
    Ybar_CT0 = mean(Y_CT0); Yvar_CT0 = var(Y_CT0)
    
    Y_CC <- rnorm(N_CC, mu_CC, sd_C)
    Ybar_CC = mean(Y_CC); Yvar_CC = var(Y_CC)
    
    data <- data.frame(N_CT = rep(N_CT,length(LOC_hetero)), N_CC = rep(N_CC,length(LOC_hetero)), 
                       N_EC1 = rep(N_EC1,length(LOC_hetero)), N_EC2 = rep(N_EC2,length(LOC_hetero)), N_EC3 = rep(N_EC3,length(LOC_hetero)),
                       Ybar_CT1 = rep(Ybar_CT1,length(LOC_hetero)), Ybar_CT0 = rep(Ybar_CT0,length(LOC_hetero)), Ybar_CC = rep(Ybar_CC,length(LOC_hetero)), 
                       Ybar_EC1 = Ybar_EC_div1, Ybar_EC2 = Ybar_EC_div2, Ybar_EC3 = Ybar_EC_div3, 
                       Yvar_CT1 = rep(Yvar_CT1,length(LOC_hetero)), Yvar_CT0 = rep(Yvar_CT0,length(LOC_hetero)), Yvar_CC = rep(Yvar_CC,length(LOC_hetero)), 
                       Yvar_EC1 = Yvar_EC_div1, Yvar_EC2 = Yvar_EC_div2, Yvar_EC3 = Yvar_EC_div3, 
                       N_CC_all = rep(N_EC1+N_EC2+N_EC3+N_CC,length(LOC_hetero)), Ybar_CC_all = rep(NA,length(LOC_hetero)), Yvar_CC_all = rep(NA,length(LOC_hetero)),
                       heterogeneity = LOC_hetero)
    
    for (n in 1:length(LOC_hetero)){
      
      Ybar_CC_all = mean(c(Y_CC,Y_EC_div1[n,],Y_EC_div2[n,],Y_EC_div3[n,])); Yvar_CC_all = var(c(Y_CC,Y_EC_div1[n,],Y_EC_div2[n,],Y_EC_div3[n,])) 
      data[n,"Ybar_CC_all"] = Ybar_CC_all; data[n,"Yvar_CC_all"] = Yvar_CC_all
    }
    
    data_norm[[i]] <- data
    
  }
  return(data_norm)
}


############################################################################################
# Conducting the simulation study of leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms
############################################################################################

######### Inputs:
#@ data: summary statistics of each simulated dataset
#@ seed: provide value for the argument of "seed" in function stan() which generates initial values for MCMC

######### Outputs:
# rst_n: the results of each simulated dataset

Norm_multiple_simulation_function <- function(data, seed){
  
  N_CT = data[1];N_CC = data[2];N_EC1 = data[3];N_EC2 = data[4];N_EC3 = data[5]; N_CC_all = data[18];
  
  Ybar_CT1 = data[6];Ybar_CT0 = data[7];Ybar_CC = data[8];Ybar_EC1 = data[9];Ybar_EC2 = data[10];Ybar_EC3 = data[11]; Ybar_CC_all = data[19];
  
  Yvar_CT1 = data[12];Yvar_CT0 = data[13];Yvar_CC = data[14];Yvar_EC1 = data[15];Yvar_EC2 = data[16];Yvar_EC3 = data[17]; Yvar_CC_all = data[20];
  
  ######### analysis
  
  fit0_n_H0    <- Norm_convention(N_CT = N_CT, Ybar_CT = Ybar_CT0, Yvar_CT = Yvar_CT0, N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, seed = seed)  %>% as.data.frame()
  fit0_n_H1    <- Norm_convention(N_CT = N_CT, Ybar_CT = Ybar_CT1, Yvar_CT = Yvar_CT1, N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, seed = seed)  %>% as.data.frame()
  
  ## Integrated Strategy
  fit1_n_I <- Norm_multi_MPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                           N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MPP1_I", seed = seed) %>% as.data.frame()
  fit2_n_I <- Norm_multi_MPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                           N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MPP2_I", seed = seed) %>% as.data.frame()
  fit3_n_I <- Norm_multi_MBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                           N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MBPP_I", seed = seed) %>% as.data.frame()
  fit4_n_I <- Norm_multi_MCBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                             N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MCBPP_I", seed = seed) %>% as.data.frame()
  
  ## Arm-based Strategy
  fit1_n_A <- Norm_multi_MPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                              N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MPP1_A", seed = seed) %>% as.data.frame()
  fit2_n_A <- Norm_multi_MPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                              N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MPP2_A", seed = seed) %>% as.data.frame()
  fit3_n_A <- Norm_multi_MBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                              N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MBPP_A", seed = seed) %>% as.data.frame()
  fit4_n_A <- Norm_multi_MCBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                                N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "MCBPP_A", seed = seed) %>% as.data.frame()
  fit5_n_A <- Norm_multi_MCBPP(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                                N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, type = "rMCBPP_A", seed = seed) %>% as.data.frame()
  
  
  fit6_n <- Norm_multi_CP_gamma(N_CC = N_CC, Ybar_CC = Ybar_CC, Yvar_CC = Yvar_CC, N_EC1 = N_EC1, Ybar_EC1 = Ybar_EC1, Yvar_EC1 = Yvar_EC1, 
                                N_EC2 = N_EC2, Ybar_EC2 = Ybar_EC2, Yvar_EC2 = Yvar_EC2, N_EC3 = N_EC3, Ybar_EC3 = Ybar_EC3, Yvar_EC3 = Yvar_EC3, seed = seed) %>% as.data.frame()
  fit7_n_H0 <- Norm_convention(N_CC = N_CC_all, Ybar_CC = Ybar_CC_all, Yvar_CC = Yvar_CC_all, N_CT = N_CT, Ybar_CT = Ybar_CT0, Yvar_CT = Yvar_CT0, seed = seed) %>% as.data.frame()
  fit7_n_H1 <- Norm_convention(N_CC = N_CC_all, Ybar_CC = Ybar_CC_all, Yvar_CC = Yvar_CC_all, N_CT = N_CT, Ybar_CT = Ybar_CT1, Yvar_CT = Yvar_CT1, seed = seed) %>% as.data.frame()
  
  
  
  ## summary-H0
  rst0_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit0_n_H0)
  
  rst1_n_H0_I <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit1_n_I)
  rst2_n_H0_I <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit2_n_I)
  rst3_n_H0_I <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit3_n_I)
  rst4_n_H0_I <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit4_n_I)
  
  rst1_n_H0_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit1_n_A)
  rst2_n_H0_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit2_n_A)
  rst3_n_H0_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit3_n_A)
  rst4_n_H0_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit4_n_A)
  rst5_n_H0_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit5_n_A)
  
  rst6_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit6_n)
  
  
  ## summary-H1
  rst0_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit0_n_H1)
  
  rst1_n_H1_I <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit1_n_I)
  rst2_n_H1_I <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit2_n_I)
  rst3_n_H1_I <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit3_n_I)
  rst4_n_H1_I <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit4_n_I)
  
  rst1_n_H1_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit1_n_A)
  rst2_n_H1_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit2_n_A)
  rst3_n_H1_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit3_n_A)
  rst4_n_H1_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit4_n_A)
  rst5_n_H1_A <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit5_n_A)
  
  rst6_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit6_n)
  
  rst7_n_H0 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H0, post_C = fit7_n_H0)
  rst7_n_H1 <- Get_normal_rst2(N_CC = N_CC, N_EC = N_EC1+N_EC2+N_EC3, post_T = fit0_n_H1, post_C = fit7_n_H1)
  
  rst_n_H0_I <- as.data.frame(rbind(rst0_n_H0,rst1_n_H0_I,rst2_n_H0_I,rst3_n_H0_I,rst4_n_H0_I,rst6_n_H0,rst7_n_H0))
  rst_n_H0_I <- data.frame(rst_n_H0_I, prior = c("No-borrowing","MPP1","MPP2","MBPP","MCBPP","CP","Pooled"), test = "H0", method = "I")
  rst_n_H1_I <- as.data.frame(rbind(rst0_n_H1,rst1_n_H1_I,rst2_n_H1_I,rst3_n_H1_I,rst4_n_H1_I,rst6_n_H1,rst7_n_H1))
  rst_n_H1_I <- data.frame(rst_n_H1_I, prior = c("No-borrowing","MPP1","MPP2","MBPP","MCBPP","CP","Pooled"), test = "H1", method = "I")
  
  rst_n_H0_A <- as.data.frame(rbind(rst0_n_H0,rst1_n_H0_A,rst2_n_H0_A,rst3_n_H0_A,rst4_n_H0_A,rst5_n_H0_A))
  rst_n_H0_A <- data.frame(rst_n_H0_A, prior = c("No-borrowing","MPP1","MPP2","MBPP","MCBPP","rMCBPP"), test = "H0", method = "A")
  rst_n_H1_A <- as.data.frame(rbind(rst0_n_H1,rst1_n_H1_A,rst2_n_H1_A,rst3_n_H1_A,rst4_n_H1_A,rst5_n_H1_A))
  rst_n_H1_A <- data.frame(rst_n_H1_A, prior = c("No-borrowing","MPP1","MPP2","MBPP","MCBPP","rMCBPP"), test = "H1", method = "A")
  
  rst_n <- rbind(rst_n_H0_I,rst_n_H1_I,rst_n_H0_A,rst_n_H1_A)
  rownames(rst_n) <- NULL
  return(rst_n)
}



