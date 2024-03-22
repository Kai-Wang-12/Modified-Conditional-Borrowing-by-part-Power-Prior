//
//  Commensurate Prior
//  To leverage external Gaussian endpoint with unknown variance from three external control arms
//  Conjugate gamma distribution as the initial prior for its random heterogeneity parameter
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
  
  //alpha
  real<lower = 0> alpha1;
  
  //beta
  real<lower = 0> beta1;

}

parameters {
  real  mu_CC; // population mean of concurrent controls
  real  mu_EC; // population variance of concurrent controls
  real<lower = 0>  sigma_CC;
  real<lower = 0>  sigma_EC1;
  real<lower = 0>  sigma_EC2;
  real<lower = 0>  sigma_EC3;
  real<lower = 0>  tau2_mu; // random heterogeneity parameter 
}

transformed parameters{
  real<lower = 0>  sigma_EC_t1;
  real<lower = 0>  sigma_EC_t2;
  real<lower = 0>  sigma_EC_t3;
  real<lower = 0>  sigma_CC_t;
  sigma_EC_t1 = (N_EC1-1)*Yvar_EC1/sigma_EC1;
  sigma_EC_t2 = (N_EC2-1)*Yvar_EC2/sigma_EC2;
  sigma_EC_t3 = (N_EC3-1)*Yvar_EC3/sigma_EC3;
  sigma_CC_t = (N_CC-1)*Yvar_CC/sigma_CC;
}

model {
  
  // Jeffreys' prior as the initial prior   
  target += -1.5*log(sigma_EC1); 
  target += -1.5*log(sigma_EC2); 
  target += -1.5*log(sigma_EC3); 
  target += -1.5*log(sigma_CC);  
  
  // Conjugate gamma distribution as the initial prior   
  tau2_mu ~ gamma(alpha1, beta1);
  
  // CP
  mu_CC ~ normal(mu_EC, sqrt(1/tau2_mu));
  
  //likelihood
  target += chi_square_lpdf(sigma_EC_t1|N_EC1-1);
  Ybar_EC1 ~ normal(mu_EC, sqrt(sigma_EC1/N_EC1));
  target += chi_square_lpdf(sigma_EC_t2|N_EC2-1);
  Ybar_EC2 ~ normal(mu_EC, sqrt(sigma_EC2/N_EC2));
  target += chi_square_lpdf(sigma_EC_t3|N_EC3-1);
  Ybar_EC3 ~ normal(mu_EC, sqrt(sigma_EC3/N_EC3));
  target += chi_square_lpdf(sigma_CC_t|N_CC-1);
  Ybar_CC ~ normal(mu_CC, sqrt(sigma_CC/N_CC));

}


