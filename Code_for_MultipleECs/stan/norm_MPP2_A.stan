//
//  Modified Power Prior
//  To leverage external Gaussian endpoint with known variance (MPP2) from three external control arms
//  Uniform distribution as the initial prior of the random discounting parameters
//  The Arm-based Strategy
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
  real<lower = 1>  N_CC1;
  real  Ybar_CC1;
  real<lower = 0>  Yvar_CC1; 
  real<lower = 1>  N_CC2;
  real  Ybar_CC2;
  real<lower = 0>  Yvar_CC2; 
  real<lower = 1>  N_CC3;
  real  Ybar_CC3;
  real<lower = 0>  Yvar_CC3;
  real<lower = 0>  sigma2_CC_hat1;
  real<lower = 0>  sigma2_CC_hat2;
  real<lower = 0>  sigma2_CC_hat3;
  real<lower = 0>  sigma2_EC_hat1;
  real<lower = 0>  sigma2_EC_hat2;
  real<lower = 0>  sigma2_EC_hat3;
  
  N_CC1 = N_CC;
  N_CC2 = N_CC;
  N_CC3 = N_CC;
  Ybar_CC1 = Ybar_CC;
  Ybar_CC2 = Ybar_CC;
  Ybar_CC3 = Ybar_CC;
  Yvar_CC1 = Yvar_CC;
  Yvar_CC2 = Yvar_CC;
  Yvar_CC3 = Yvar_CC;
  sigma2_EC_hat1 = (N_EC1-1)*Yvar_EC1/N_EC1;
  sigma2_EC_hat2 = (N_EC2-1)*Yvar_EC2/N_EC2;
  sigma2_EC_hat3 = (N_EC3-1)*Yvar_EC3/N_EC3;
  sigma2_CC_hat1 = (N_CC1-1)*Yvar_CC1/N_CC1;
  sigma2_CC_hat2 = (N_CC2-1)*Yvar_CC2/N_CC2;
  sigma2_CC_hat3 = (N_CC3-1)*Yvar_CC3/N_CC3;
  
}

parameters {
  real<lower = 0, upper = 1>  delta1;  // random power parameter specific to EC1
  real<lower = 0, upper = 1>  delta2;  // random power parameter specific to EC2
  real<lower = 0, upper = 1>  delta3;  // random power parameter specific to EC3
  real  mu_CC1;                                         // population mean of concurrent controls specific to EC1
  real  mu_CC2;                                         // population mean of concurrent controls specific to EC2
  real  mu_CC3;                                         // population mean of concurrent controls specific to EC3
}

model {
  
  // Uniform distribution as the initial prior
  delta1 ~ beta(1, 1);
  delta2 ~ beta(1, 1);
  delta3 ~ beta(1, 1);
  
  // MPP2 constructed using each external control arm and the Jeffreys' prior
  target += normal_lpdf(mu_CC1|Ybar_EC1, sqrt(sigma2_EC_hat1/(N_EC1*delta1)));
  target += normal_lpdf(mu_CC2|Ybar_EC2, sqrt(sigma2_EC_hat2/(N_EC2*delta2)));
  target += normal_lpdf(mu_CC3|Ybar_EC3, sqrt(sigma2_EC_hat3/(N_EC3*delta3)));
  
  //likelihood
  Ybar_CC1 ~ normal(mu_CC1, sqrt(sigma2_CC_hat1/N_CC1));
  Ybar_CC2 ~ normal(mu_CC2, sqrt(sigma2_CC_hat2/N_CC2));
  Ybar_CC3 ~ normal(mu_CC3, sqrt(sigma2_CC_hat3/N_CC3));
  
}

generated quantities{
  real<lower = 0>  sigma2_CC_hat;
  real<lower = 0>  sigma_all;
  real mu_CC;     // population mean of concurrent controls
  
  // posterior distribution
  sigma2_CC_hat = (N_CC-1)*Yvar_CC/N_CC;
  sigma_all = 1/(delta1*N_EC1/sigma2_EC_hat1 + delta2*N_EC2/sigma2_EC_hat2 + delta3*N_EC3/sigma2_EC_hat3 + N_CC/sigma2_CC_hat);
  mu_CC = normal_rng((delta1*N_EC1*Ybar_EC1/sigma2_EC_hat1 + delta2*N_EC2*Ybar_EC2/sigma2_EC_hat2 + delta3*N_EC3*Ybar_EC3/sigma2_EC_hat3 + N_CC*Ybar_CC/sigma2_CC_hat)*sigma_all, sqrt(sigma_all));
}

