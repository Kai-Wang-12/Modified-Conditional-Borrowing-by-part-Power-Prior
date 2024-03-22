//
//  Modified Power Prior
//  To leverage external Gaussian endpoint with known variance (MPP2) from three external control arms
//  Uniform distribution as the initial prior of the random discounting parameters
//  The Integrated Strategy
//

data {
  
  //external controls (EC1, EC2, EC3)
  real<lower = 1>  N_EC1;
  real  Ybar_EC1;
  real<lower = 0>  Yvar_EC1; 
  
  real<lower = 1>  N_EC2;
  real  Ybar_EC2;
  real<lower = 0>  Yvar_EC2; 
  
  real<lower = 1>  N_EC3;
  real  Ybar_EC3;
  real<lower = 0>  Yvar_EC3; 
  
  //concurrent controls
  real<lower = 1>  N_CC;
  real  Ybar_CC;
  real<lower = 0>  Yvar_CC; 

}

transformed data{
  real<lower = 0>  sigma2_EC_hat1;
  real<lower = 0>  sigma2_EC_hat2;
  real<lower = 0>  sigma2_EC_hat3;
  real<lower = 0>  sigma2_CC_hat;
  sigma2_EC_hat1 = (N_EC1-1)*Yvar_EC1/N_EC1;
  sigma2_EC_hat2 = (N_EC2-1)*Yvar_EC2/N_EC2;
  sigma2_EC_hat3 = (N_EC3-1)*Yvar_EC3/N_EC3;
  sigma2_CC_hat = (N_CC-1)*Yvar_CC/N_CC;
  
}

parameters {
  real<lower = 0, upper = 1>  delta1;  // random power parameter specific to EC1
  real<lower = 0, upper = 1>  delta2;  // random power parameter specific to EC2
  real<lower = 0, upper = 1>  delta3;  // random power parameter specific to EC3
  real  mu_CC;                                            // population mean of concurrent controls
}

transformed parameters{
  real<lower = 0>  sigma_EC_all;
  sigma_EC_all = 1/(delta1*N_EC1/sigma2_EC_hat1 + delta2*N_EC2/sigma2_EC_hat2 + delta3*N_EC3/sigma2_EC_hat3);
}

model {
  
  // Uniform distribution as the initial prior
  delta1 ~ beta(1, 1);
  delta2 ~ beta(1, 1);
  delta3 ~ beta(1, 1);
  
  // MPP2 constructed using all external controls and the Jeffreys' prior
  target += normal_lpdf(mu_CC|(delta1*N_EC1*Ybar_EC1/sigma2_EC_hat1 + delta2*N_EC2*Ybar_EC2/sigma2_EC_hat2 + delta3*N_EC3*Ybar_EC3/sigma2_EC_hat3)*sigma_EC_all, sqrt(sigma_EC_all));
  
  //likelihood
  Ybar_CC ~ normal(mu_CC, sqrt(sigma2_CC_hat/N_CC));
  
}

