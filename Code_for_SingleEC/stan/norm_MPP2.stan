//
//  Modified Power Prior
//  To leverage external Gaussian endpoint with known variance (MPP2) from a single external control arm
//  Uniform distribution as the initial prior of the single random discounting parameter
//

data {
  
  //external controls
  real<lower = 1>  N_EC;
  real  Ybar_EC;
  real<lower = 0>  Yvar_EC; 
  
  //concurrent controls
  real<lower = 1>  N_CC;
  real  Ybar_CC;
  real<lower = 0>  Yvar_CC; 

}

transformed data{
  real<lower = 0>  sigma2_EC_hat;
  real<lower = 0>  sigma2_CC_hat;
  sigma2_EC_hat = (N_EC-1)*Yvar_EC/N_EC;
  sigma2_CC_hat = (N_CC-1)*Yvar_CC/N_CC;
  
}

parameters {
  real<lower = 0, upper = 1>  delta;    // random power parameter
  real  mu_CC;                                            // population mean of concurrent controls
}

model {
  
  // Uniform distribution as the initial prior
  delta ~ beta(1,1);
  
  // MPP2 constructed using external controls and the Jeffreys' prior
  target += normal_lpdf(mu_CC|Ybar_EC, sqrt(sigma2_EC_hat/(N_EC*delta)));
  
  //likelihood
  Ybar_CC ~ normal(mu_CC, sqrt(sigma2_CC_hat/N_CC));
  
}

