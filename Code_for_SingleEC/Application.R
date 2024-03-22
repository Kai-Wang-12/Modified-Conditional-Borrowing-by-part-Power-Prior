##################################################################################################################################################################
#-------------------------------------------       This file includes R codes to conduct the application      ---------------------------------------------------#
##################################################################################################################################################################

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

############################################################################################
#-----------------------------  Simulation_for_Calibration  -------------------------------#
############################################################################################

#----------------------------------- Warm-up running  

sd_C = 2.5

fit0_n <- Norm_convention(N_CT = 202, Ybar_CT = 0, Yvar_CT = sd_C^2, N_CC = 101, Ybar_CC = 0, Yvar_CC = sd_C^2, seed = seed)
fit1_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MPP1", seed = seed)
fit2_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MPP2", seed = seed)
fit3_n <- Norm_MBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MBPP", seed = seed)
fit4_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MCBPP", seed = seed)
fit5_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "rMCBPP", seed = seed)
fit6_n <- Norm_CP_gamma(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, seed = seed)

#----------------------------------- dataset 

Creat_dataset_normal_application(N_CC = 139, N_CT = 280, N_EC = 271, 
                                 norm_effect = 0, mu_EC = -0.87, mu_CC = -0.87,  # under H0 and none location heterogeneity (mu_CC = mu_EC)
                                 sd_CT = 3.333, sd_CC = 2.833, sd_EC = 2.833,
                                 N_dataset = 2500, seed=202306) %>%
  saveRDS(paste0(wd,"/simulation_dataset/data_norm_E1.rds"))

Creat_dataset_normal_application(N_CC = 139, N_CT = 280, N_EC = 271,
                                 norm_effect = 0, mu_EC = -0.87, mu_CC = -0.87,  # under H0 and none location heterogeneity (mu_CC = mu_EC)
                                 sd_CT = 3.333, sd_CC = 2.833, sd_EC = 3.004,
                                 N_dataset = 2500, seed=202306) %>%
  saveRDS(paste0(wd,"/simulation_dataset/data_norm_E2.rds"))

Creat_dataset_normal_application(N_CC = 139, N_CT = 280, N_EC = 271,
                                 norm_effect = 0, mu_EC = -0.87, mu_CC = -0.87,  # under H0 and none location heterogeneity (mu_CC = mu_EC)
                                 sd_CT = 3.333, sd_CC = 2.833, sd_EC = 5,
                                 N_dataset = 2500, seed=202306) %>%
  saveRDS(paste0(wd,"/simulation_dataset/data_norm_E3.rds"))


#----------------------------------- Simulation 

data_name <- paste0("data_norm_E",c(1:3),".rds")
norm_effect = 0

for(i in 1:length(data_name)){
  
  dataset_name <- data_name[i]
  dataset <- readRDS(paste0(wd,"/simulation_dataset/",dataset_name))
  
  apply(dataset, 1, function(y){
    Norm_simulation_function(as.matrix(y), Pooled = F, seed = seed) %>%   
      data.table()
  })%>%
    rbindlist() %>%
    data.table() %>%
    saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_E",i,".rds"))
}


for(i in 1:length(data_name)){
  
  data <- readRDS((file = paste0(wd, "/rst_for_SingleEC/simu_E",i,".rds")))
  
  threshold <- filter(data, test=="H0") %>%                                         # Find threshold value under H0 and none location heterogeneity
    group_by(prior) %>%
    summarise(threshold = quantile(1 - post_probability, probs = 0.975)) %>%        # H1: ATE<0;  Control type I error rate at 0.025
    saveRDS(file = paste0(wd, "/rst_for_SingleEC/simu_E_threshold_",i,".rds"))
  }
gc()

rst1_H0 <- read_rds(file = paste0(wd, "/rst_for_SingleEC/simu_E_threshold_1.rds"))
rst2_H0 <- read_rds(file = paste0(wd, "/rst_for_SingleEC/simu_E_threshold_2.rds"))
rst3_H0 <- read_rds(file = paste0(wd, "/rst_for_SingleEC/simu_E_threshold_3.rds"))


############################################################################################
#-------------------------------------- Analysis  -----------------------------------------#
############################################################################################

data1 <- c(280, 139, 271, -2.65, -0.87, -1.03, 3.333^2, 2.833^2, 2.833^2)
rst_1 <- Norm_application_function(data1, seed = 202306, norm_effect = -1.78)
rst_1

data2 <- c(280, 139, 271, -2.65, -0.87, -1.03, 3.333^2, 2.833^2, 3.004^2)
rst_2 <- Norm_application_function(data2, seed = 202306, norm_effect = -1.78)
rst_2

data3 <- c(280, 139, 271, -2.65, -0.87, -1.03, 3.333^2, 2.833^2, 5^2)
rst_3 <- Norm_application_function(data3, seed = 202306, norm_effect = -1.78)
rst_3

norm_effect = -1.78

rst_case1 <- rst_1 %>%
  full_join(rst1_H0, by = "prior") %>%
  group_by(prior) %>%
  summarise(ESS = ESS, Proportion = Proportion, Est = Delta_mean,
            Bias = Delta_mean - norm_effect, 
            MSE = MSE, 
            post_prob = 1 - post_probability,               # H1: ATE<0
            threshold = threshold,
            Success = post_prob > 0.975,                    # Success based on 0.975
            Success_calibrated = post_prob > threshold,     # calibrated Success based on threshold
            CI_length = CI_length) %>%
  as.data.frame()  %>%
  saveRDS(file = paste0(wd,"/rst_for_SingleEC/rst_case1.rds"))

rst_case2 <- rst_2 %>%
  full_join(rst2_H0, by = "prior") %>%
  group_by(prior) %>%
  summarise(ESS = ESS, Proportion = Proportion, Est = Delta_mean,
            Bias = Delta_mean - norm_effect, 
            MSE = MSE, 
            post_prob = 1 - post_probability,               # H1: ATE<0 
            threshold = threshold,
            Success = post_prob > 0.975,                    # Success based on 0.975
            Success_calibrated = post_prob > threshold,     # calibrated Success based on threshold
            CI_length = CI_length) %>%
  as.data.frame() %>%
  saveRDS(file = paste0(wd,"/rst_for_SingleEC/rst_case2.rds"))

rst_case3 <- rst_3 %>%
  full_join(rst3_H0, by = "prior") %>%
  group_by(prior) %>%
  summarise(ESS = ESS, Proportion = Proportion, Est = Delta_mean,
            Bias = Delta_mean - norm_effect, 
            MSE = MSE, 
            post_prob = 1 - post_probability,               # H1: ATE<0
            threshold = threshold,
            Success = post_prob > 0.975,                    # Success based on 0.975
            Success_calibrated = post_prob > threshold,     # calibrated Success based on threshold
            CI_length = CI_length) %>%
  as.data.frame() %>%
  saveRDS(file = paste0(wd,"/rst_for_SingleEC/rst_case3.rds"))



