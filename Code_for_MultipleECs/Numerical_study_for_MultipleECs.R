#########################################################################################################################################################################################################
#------------     This file includes R codes to conduct the numerical study of leveraging external Gaussian endpoint with unknown variance from multiple (three) external control arms      ------------#
#########################################################################################################################################################################################################

rm(list=ls())

############################################################################################
#---------------------------- Loading packages and functions  -----------------------------#
############################################################################################

library(tidyverse)
library(dplyr)
library(data.table)
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
VAR_hetero = c(seq(0.1,2,0.1),seq(2.25,5,0.25))                                   # Variance heterogeneity
sd_C <- 2.5;

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
#----------------------------- Numerical study for Figure 2  ------------------------------#
############################################################################################

#------------- Situation1  ---------------#
lapply(VAR_hetero, function(x){
  Norm_multiple_numerical_function(N = c(426, 213, 71, 71, 71), Ybar = c(0, 0, 0, 1*sd_C, -1*sd_C), Yvar = c(sd_C^2, sd_C^2, sd_C^2, sd_C^2 * (x)^2, sd_C^2 * (x)^2), seed = seed) %>%
    data.frame(VAR_hetero = x)
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/rst_for_MultipleECs/num_A.rds"))

#------------- Situation2  ---------------#
lapply(VAR_hetero, function(x){
  Norm_multiple_numerical_function(N = c(426, 213, 71, 71, 71), Ybar = c(0, 0, 0, 0, 0), Yvar = c(sd_C^2, sd_C^2, sd_C^2, sd_C^2 * (x)^2, sd_C^2 * (x)^2), seed = seed) %>%
    data.frame(VAR_hetero = x)
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/rst_for_MultipleECs/num_B.rds"))

gc()

