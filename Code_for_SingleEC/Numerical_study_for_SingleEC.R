#########################################################################################################################################################################################################
#-------------        This file includes R codes to conduct the numerical study of leveraging external Gaussian endpoint with unknown variance from a single external control arm        ---------------#
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

setwd("~/Code_for_SingleEC")  #  set the working directory
source("Function_for_SingleEC.R")

############################################################################################
#----------------------------------- General setting  ------------------------------------#
############################################################################################
wd <- getwd()
setwd("stan")  # Set the working directory to the stan code file
seed =2023
LOC_hetero = seq(-2,2,0.1)                                                      # Location heterogeneity
VAR_hetero = c(seq(0.1,2,0.1),seq(2.25,5,0.25))                                 # Variance heterogeneity
sd_CC <- 2.5;  # sample standard error of concurrent controls

############################################################################################
#----------------------------------- Warm-up running  -------------------------------------#
############################################################################################

fit0_n <- Norm_convention(N_CT = 202, Ybar_CT = 0, Yvar_CT = sd_CC^2, N_CC = 101, Ybar_CC = 0, Yvar_CC = sd_CC^2, seed = seed)
fit1_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "MPP1", seed = seed)
fit2_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "MPP2", seed = seed)
fit3_n <- Norm_MBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "MBPP", seed = seed)
fit4_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "MCBPP", seed = seed)
fit5_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "rMCBPP", seed = seed)
fit6_n <- Norm_CP_gamma(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, seed = seed)

gc()

############################################################################################
#-------------------------- Numerical study for Figure 1, S2-4  ---------------------------#
############################################################################################

N = c(202, 101, 101)   # Figure1
#N = c(26, 13, 13)      # FigureS2
#N = c(202, 152, 50)    # FigureS3
#N = c(202, 50, 152)    # FigureS4

#------------- None location heterogeneity  ---------------#
lapply(VAR_hetero, function(x){
  Norm_numerical_function(N = c(202, 101, 101), Ybar = c(0, 0, 0), Yvar = c(sd_CC^2, sd_CC^2, sd_CC^2 * (x)^2), seed = seed) %>%
    data.frame(VAR_hetero = x)
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/rst_for_SingleEC/num_A1.rds"))

#------------- None variance heterogeneity  ---------------#
lapply(LOC_hetero, function(x){
  Norm_numerical_function(N = c(202, 101, 101), Ybar = c(0, 0, x*sd_CC), Yvar = c(sd_CC^2, sd_CC^2, sd_CC^2), seed = seed) %>%
    data.frame(LOC_hetero = x) 
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/rst_for_SingleEC/num_A4.rds"))

#------------- Mild variance heterogeneity  ---------------#
lapply(LOC_hetero, function(x){
  Norm_numerical_function(N = c(202, 101, 101), Ybar = c(0, 0, x*sd_CC), Yvar = c(sd_CC^2, sd_CC^2, 1.2^2*sd_CC^2), seed = seed) %>%
    data.frame(LOC_hetero = x) 
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/rst_for_SingleEC/num_A5.rds"))
lapply(LOC_hetero, function(x){
  Norm_numerical_function(N = c(202, 101, 101), Ybar = c(0, 0, x*sd_CC), Yvar = c(sd_CC^2, sd_CC^2, 0.8^2*sd_CC^2), seed = seed) %>%
    data.frame(LOC_hetero = x) 
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/rst_for_SingleEC/num_A6.rds"))

#------------- Severe variance heterogeneity  ---------------#
lapply(LOC_hetero, function(x){
  Norm_numerical_function(N = c(202, 101, 101), Ybar = c(0, 0, x*sd_CC), Yvar = c(sd_CC^2, sd_CC^2, 2.5^2*sd_CC^2), seed = seed) %>%
    data.frame(LOC_hetero = x) 
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/rst_for_SingleEC/num_A7.rds"))
lapply(LOC_hetero, function(x){
  Norm_numerical_function(N = c(202, 101, 101), Ybar = c(0, 0, x*sd_CC), Yvar = c(sd_CC^2, sd_CC^2, 0.1^2*sd_CC^2), seed = seed) %>%
    data.frame(LOC_hetero = x) 
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/rst_for_SingleEC/num_A8.rds"))

gc()
