#########################################################################################################################################################################################################
#-----------     This file includes R codes to conduct the simulation study of leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms      ------------#
#########################################################################################################################################################################################################

rm(list=ls())

############################################################################################
#---------------------------- Loading packages and functions  -----------------------------#
############################################################################################

library(tidyverse)
library(data.table)
library(dplyr)
library(rlist)
library(rstan)
library(rstantools)
library(RBesT)
library(MASS)
library(corpcor)
library(Hmisc)

setwd("~/Code_for_MultipleECs")  #  set the working directory
source("Function_for_MultipleECs.R")

############################################################################################
#----------------------------------- General setting  ------------------------------------#
############################################################################################
wd <- getwd()
setwd("stan")  # Set the working directory to the stan code file
seed =2023
LOC_hetero = seq(-2,2,0.1)                                                # Location heterogeneity
VAR_hetero = c(1,0.1,2.5)                                                 # Variance heterogeneity
sd_C = 2.5

############################################################################################
#---------- Generating simulated datasets for the multiple external control arms-----------#
############################################################################################

norm_effect = 0.48

for (i in 1:3){
  Creat_dataset_MultipleECs(N_CC = 213, N_CT = 426, N_EC1 = 71, N_EC2 = 71, N_EC3 = 71, 
                            mu_CC = 0, mu_EC1 = 0, mu_EC2 = sd_C, mu_EC3 = -1*sd_C,
                            LOC_hetero = LOC_hetero, norm_effect = norm_effect, 
                            sd_C = sd_C, sd_EC1 = sd_C, sd_EC2 = VAR_hetero[i]*sd_C, sd_EC3 = VAR_hetero[i]*sd_C, 
                            N_dataset = 2500, seed = 202305) %>%
    rbindlist() %>%
    saveRDS(paste0(wd,"/simulation_dataset/data_norm_A",i,".rds")) ##  Scenario I-K
}

for (i in 1:3){
  Creat_dataset_MultipleECs(N_CC = 213, N_CT = 426, N_EC1 = 71, N_EC2 = 71, N_EC3 = 71, 
                            mu_CC = 0, mu_EC1 = 0, mu_EC2 = 0, mu_EC3 = 0,
                            LOC_hetero = LOC_hetero, norm_effect = norm_effect, 
                            sd_C = sd_C, sd_EC1 = sd_C, sd_EC2 = VAR_hetero[i]*sd_C, sd_EC3 = VAR_hetero[i]*sd_C, 
                            N_dataset = 2500, seed = 202305) %>%
    rbindlist() %>%
    saveRDS(paste0(wd,"/simulation_dataset/data_norm_B",i,".rds")) ##  Scenario F-H
}


############################################################################################
#----------------------------------- Warm-up running  -------------------------------------#
############################################################################################
fit0_n <- Norm_convention(N_CT = 426, Ybar_CT = 0, Yvar_CT = sd_C^2, N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, seed = seed)


fit1_n_I <- Norm_multi_MPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                           N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MPP1_I", seed = seed)
fit1_n_A <- Norm_multi_MPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                           N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MPP1_A", seed = seed) 
fit2_n_I <- Norm_multi_MPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                           N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MPP2_I", seed = seed) 
fit2_n_A <- Norm_multi_MPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                           N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MPP2_A", seed = seed)
fit3_n_I <- Norm_multi_MBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                            N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MBPP_I", seed = seed) 
fit3_n_A <- Norm_multi_MBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                            N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MBPP_A", seed = seed) 
fit4_n_I <- Norm_multi_MCBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                             N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MCBPP_I", seed = seed) 
fit4_n_A <- Norm_multi_MCBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                             N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MCBPP_A", seed = seed) 
fit5_n_A <- Norm_multi_MCBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                             N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "rMCBPP_A", seed = seed) 
fit6_n <- Norm_multi_CP_gamma(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                              N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, seed = seed)

gc()

############################################################################################
#----------------------------------- Simulation study -------------------------------------#
############################################################################################

data_name <- paste0("data_norm_A",c(1:3),".rds");  norm_effect = 0.48  ##  Scenario I-K
# data_name <- paste0("data_norm_B",c(1:3),".rds");  norm_effect = 0.48  ##  Scenario F-H

for(i in 1:length(data_name)){
  
  apply(readRDS(paste0(wd,"/simulation_dataset/",data_name[i])), 1, function(y){
    Norm_multiple_simulation_function(y, seed = seed) %>%   
      data.table(heterogeneity = y[21])
  })%>%
    rbindlist() %>%
    data.table() %>%
    saveRDS(file = paste0(wd, "/rst_for_MultipleECs/simu_A",i,".rds")) ##  Scenario I-K
  # saveRDS(file = paste0(wd, "/rst_for_MultipleECs/simu_B",i,".rds")) ##  Scenario F-H
}

norm_effect = 0.48

for(i in 1:length(data_name)){
  
  threshold <- filter(readRDS(file = paste0(wd, "/rst_for_MultipleECs/simu_B",i,".rds")), 
                      test=="H0", heterogeneity == 0) %>%                       # Find threshold value under H0 and none location heterogeneity (Scenario F-H)
    group_by(prior,method) %>%
    summarise(threshold = quantile(post_probability, probs = 0.975))            # Control type I error rate at 0.025
  
  data <- readRDS(file = paste0(wd, "/rst_for_MultipleECs/simu_A",i,".rds"))
  #data <- readRDS((file = paste0(wd, "/rst_for_MultipleECs/simu_B",i,".rds")))
  
  rst_H0 <- data %>%
    filter(test=="H0") %>%   
    full_join(threshold, by = c("prior","method")) %>%
    group_by(prior,method,heterogeneity) %>%
    summarise(ESS = median(ESS), Proportion = median(Proportion), 
              Bias = mean(Delta_mean - 0), 
              MSE = mean((Delta_mean - 0)^2), 
              Type_I_error = mean(post_probability > 0.975),                    # Type I error rate based on 0.975
              Type_I_error_calibrated = mean(post_probability > threshold),     # calibrated Type I error rate based on threshold
              CI_length = mean(CI_length), test = "H0") %>%
    arrange(prior, method, heterogeneity) %>%
    as.data.frame() %>%
    saveRDS(file = paste0(wd, "/rst_for_MultipleECs/simu_A_rst_H0_",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_MultipleECs/simu_B_rst_H0_",i,".rds"))
  
  rst_H1 <- data %>%
    filter(test=="H1") %>%   
    full_join(threshold, by = c("prior","method")) %>%
    group_by(prior,method,heterogeneity) %>%
    summarise(ESS = median(ESS), Proportion = median(Proportion), 
              Bias = mean(Delta_mean - norm_effect), 
              MSE = mean((Delta_mean - norm_effect)^2), 
              Power = mean(post_probability > 0.975),                           # Power based on 0.975
              Power_calibrated = mean(post_probability > threshold),            # calibrated Power based on threshold
              CI_length = mean(CI_length), test = "H1") %>%
    arrange(prior, method, heterogeneity) %>%
    as.data.frame() %>%
    saveRDS(file = paste0(wd, "/rst_for_MultipleECs/simu_A_rst_H1_",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_MultipleECs/simu_B_rst_H1_",i,".rds"))
  
}
gc()



