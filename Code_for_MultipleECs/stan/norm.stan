//
//  Conventional Bayesian analysis
//  To leverage external Gaussian endpoint with unknown variance from three external control arms
//

data {
  
  //external controls (EC1, EC2, EC3)
  int<lower = 1> N_CC;
  real Ybar_CC;
  real<lower = 0> Yvar_CC;
  
  //experiment arm
  real<lower = 1> N_CT;
  real Ybar_CT;
  real<lower = 0> Yvar_CT;
}

parameters {
  real  mu_CC;                             // population mean of concurrent controls
  real<lower = 0>  sigma_CC; // population variance of concurrent controls
  real  mu_CT;                             // population mean of experiment arm
  real<lower = 0>  sigma_CT; // population variance of experiment arm
}

transformed parameters {
  real<lower = 0>  sigma_CT_t;
  real<lower = 0>  sigma_CC_t;
  sigma_CT_t = (N_CT-1)*Yvar_CT/sigma_CT;
  sigma_CC_t = (N_CC-1)*Yvar_CC/sigma_CC;
}

model {
  
  //  Jeffreys' prior as the initial prior
  target += - 1.5 * log(sigma_CT);  
  target += - 1.5 * log(sigma_CC);  
  
  //  likelihood
  target += chi_square_lpdf(sigma_CT_t|N_CT-1);
  Ybar_CT ~ normal(mu_CT, sqrt(sigma_CT/N_CT));
  
  target += chi_square_lpdf(sigma_CC_t|N_CC-1);
  Ybar_CC ~ normal(mu_CC, sqrt(sigma_CC/N_CC));
}

