//
//  Modified Conditional Borrowing-by-part Power Prior (MCBPP)
//  To leverage external Gaussian endpoint with unknown variance from a single external control arm
//  Uniform distribution as the initial prior of the random discounting parameters
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

parameters {
  real<lower = 0, upper = 1>  delta_mu;      // random mean-discounting parameter
  real<lower = 0, upper = 1>  delta_sigma; // random variance-discounting parameter
  real  mu_CC;                                                     // population mean of concurrent controls
  real<lower = 0>  sigma_CC;                         // population variance of concurrent controls
}

transformed parameters{
  real<lower = 0>  sigma_CC_t;
  sigma_CC_t = (N_CC-1)*Yvar_CC/sigma_CC;
}

model {
  
  // Uniform distribution as the initial prior
  delta_mu ~ beta(1, 1);
  delta_sigma ~ beta(1, 1);
  
  // MCBPP constructed using external controls and the Jeffreys' prior
  sigma_CC ~ inv_gamma((delta_sigma*(N_EC-1)+1)/2, delta_sigma*(N_EC-1)*Yvar_EC/2);
  target += normal_lpdf(mu_CC|Ybar_EC, sqrt(sigma_CC/(N_EC*delta_mu)));
  
  //likelihood
  target += chi_square_lpdf(sigma_CC_t|N_CC-1);
  Ybar_CC ~ normal(mu_CC, sqrt(sigma_CC/N_CC));
  
}

