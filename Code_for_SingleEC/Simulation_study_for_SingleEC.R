#########################################################################################################################################################################################################
#--------------   This file includes R codes to conduct the simulation study of leveraging external Gaussian endpoint with unknown variance from a single external control arm     ---------------------#
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

setwd("~/Code_for_SingleEC")  #  set the working directory
source("Function_for_SingleEC.R")

############################################################################################
#----------------------------------- General setting  ------------------------------------#
############################################################################################
wd <- getwd()
setwd("stan")  # Set the working directory to the stan code file
seed =2023
LOC_hetero = seq(-2,2,0.1)            # Location heterogeneity
VAR_hetero = c(1,1.2,0.8,2.5,0.1)     # Variance heterogeneity


############################################################################################
#----------- Generating simulated datasets for the single external control arm-------------#
############################################################################################

mu_CC = 0; sd_C = 2.5

N_CT = 202; N_CC = 101; N_EC = 101; norm_effect = 0.7   ##  Design 1
#N_CT = 26; N_CC = 13;  N_EC = 13;  norm_effect = 2.0   ##  Design 2
#N_CT = 202; N_CC = 152; N_EC = 50; norm_effect = 0.7   ##  Design 3
#N_CT = 202; N_CC = 50; N_EC = 152; norm_effect = 0.7   ##  Design 4

for (i in 1:5){
  Creat_dataset_SingleEC(N_CC = N_CC, N_CT = N_CT, N_EC = N_EC, 
                         norm_effect = norm_effect, mu_CC = mu_CC, sd_C = sd_C, 
                         LOC_hetero = LOC_hetero, VAR_hetero = VAR_hetero[i],
                         N_dataset = 2500, seed=202301) %>%
    rbindlist() %>%
    saveRDS(paste0(wd,"/simulation_dataset/data_norm_A",i,".rds"))
  # saveRDS(paste0(wd,"/simulation_dataset/data_norm_B",i,".rds")) ##  Design 2, seed=202302
  # saveRDS(paste0(wd,"/simulation_dataset/data_norm_C",i,".rds")) ##  Design 3, seed=202303
  # saveRDS(paste0(wd,"/simulation_dataset/data_norm_D",i,".rds")) ##  Design 4, seed=202304
}


############################################################################################
#----------------------------------- Warm-up running  -------------------------------------#
############################################################################################

fit0_n <- Norm_convention(N_CT = 202, Ybar_CT = 0, Yvar_CT = sd_C^2, N_CC = 101, Ybar_CC = 0, Yvar_CC = sd_C^2, seed = seed)
fit1_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MPP1", seed = seed)
fit2_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MPP2", seed = seed)
fit3_n <- Norm_MBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MBPP", seed = seed)
fit4_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MCBPP", seed = seed)
fit5_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "rMCBPP", seed = seed)
fit6_n <- Norm_CP_gamma(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, seed = seed)


############################################################################################
#----------------------------------- Simulation study -------------------------------------#
############################################################################################

data_name <- paste0("data_norm_A",c(1:5),".rds");  norm_effect = 0.7 ##  Design 1
# data_name <- paste0("data_norm_B",c(1:5),".rds");  norm_effect = 2.0 ##  Design 2
# data_name <- paste0("data_norm_C",c(1:5),".rds");  norm_effect = 0.7 ##  Design 3
# data_name <- paste0("data_norm_D",c(1:5),".rds");  norm_effect = 0.7 ##  Design 4

for(i in 1:length(data_name)){
  
  apply(readRDS(paste0(wd,"/simulation_dataset/",data_name[i])), 1, function(y){
    Norm_simulation_function(y, seed = seed) %>%   
      data.table(heterogeneity = y[15])
  })%>%
    rbindlist() %>%
    data.table() %>%
    saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_A",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_B",i,".rds")) ##  Design 2
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_C",i,".rds")) ##  Design 3
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_D",i,".rds")) ##  Design 4
}


for(i in 1:length(data_name)){
  
  data <- readRDS((file = paste0(wd, "/rst_for_SingleEC/simu_A",i,".rds")))
  # data <- readRDS((file = paste0(wd, "/rst_for_SingleEC/simu_B",i,".rds")))
  # data <- readRDS((file = paste0(wd, "/rst_for_SingleEC/simu_C",i,".rds")))
  # data <- readRDS((file = paste0(wd, "/rst_for_SingleEC/simu_D",i,".rds")))
  
  threshold <- filter(data, test=="H0", heterogeneity==0) %>%                   # Find threshold value under H0 and none location heterogeneity
    group_by(prior) %>%
    summarise(threshold = quantile(post_probability, probs = 0.975))            # Control type I error rate at 0.025
  
  rst_H0 <- data %>%
    filter(test=="H0") %>%   
    full_join(threshold, by = "prior") %>%
    group_by(prior,heterogeneity) %>%
    summarise(ESS = median(ESS), Proportion = median(Proportion), 
              Bias = mean(Delta_mean - 0), 
              MSE = mean((Delta_mean - 0)^2), 
              Type_I_error = mean(post_probability > 0.975),                    # Type I error rate based on 0.975
              Type_I_error_calibrated = mean(post_probability > threshold),     # calibrated Type I error rate based on threshold
              CI_length = mean(CI_length), test = "H0") %>%
    arrange(prior, heterogeneity) %>%
    as.data.frame() %>%
    saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_A_rst_H0_",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_B_rst_H0_",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_C_rst_H0_",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_D_rst_H0_",i,".rds"))
  
  rst_H1 <- data %>%
    filter(test=="H1") %>%   
    full_join(threshold, by = "prior") %>%
    group_by(prior,heterogeneity) %>%
    summarise(ESS = median(ESS), Proportion = median(Proportion), 
              Bias = mean(Delta_mean - norm_effect), 
              MSE = mean((Delta_mean - norm_effect)^2), 
              Power = mean(post_probability > 0.975),                           # Power based on 0.975
              Power_calibrated = mean(post_probability > threshold),            # calibrated Power based on threshold
              CI_length = mean(CI_length), test = "H1") %>%
    arrange(prior, heterogeneity) %>%
    as.data.frame() %>%
    saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_A_rst_H1_",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_B_rst_H1_",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_C_rst_H1_",i,".rds"))
  # saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_D_rst_H1_",i,".rds"))
  
}
gc()



