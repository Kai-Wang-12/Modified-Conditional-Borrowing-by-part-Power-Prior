//
//  Commensurate Prior
//  To leverage external Gaussian endpoint with unknown variance from a single external control arm
//  Conjugate gamma distribution as the initial prior for its random heterogeneity parameter
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

  //alpha
  real<lower = 0> alpha1;
  
  //beta
  real<lower = 0> beta1;
  
}

parameters {
  real  mu_CC; // population mean of concurrent controls
  real  mu_EC; // population mean of external controls
  real<lower = 0>  sigma_CC; // population variance of concurrent controls
  real<lower = 0>  sigma_EC; // population variance of external controls
  real<lower = 0>  tau2_mu;  // random heterogeneity parameter 
}

transformed parameters{
  real<lower = 0>  sigma_EC_t;
  real<lower = 0>  sigma_CC_t;
  sigma_EC_t = (N_EC-1)*Yvar_EC/sigma_EC;
  sigma_CC_t = (N_CC-1)*Yvar_CC/sigma_CC;
}

model {
  
  // Jeffreys' prior as the initial prior   
  target += -1.5 * log(sigma_EC);
  target += -1.5 * log(sigma_CC);
  
  // Conjugate gamma distribution as the initial prior
  tau2_mu ~ gamma(alpha1, beta1);
  
  // CP
  mu_CC ~ normal(mu_EC, sqrt(1/tau2_mu));
  
  //likelihood
  target += chi_square_lpdf(sigma_EC_t|N_EC-1);
  Ybar_EC ~ normal(mu_EC, sqrt(sigma_EC/N_EC));
  target += chi_square_lpdf(sigma_CC_t|N_CC-1);
  Ybar_CC ~ normal(mu_CC, sqrt(sigma_CC/N_CC));

}

