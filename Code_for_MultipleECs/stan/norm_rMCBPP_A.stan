//
//  Robust Modified Conditional Borrowing-by-part Power Prior (rMCBPP) with all variance-discounting parameters fixed at zero
//  To leverage external Gaussian endpoint with unknown variance from three external control arms
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
  
  N_CC1 = N_CC;
  N_CC2 = N_CC;
  N_CC3 = N_CC;
  Ybar_CC1 = Ybar_CC;
  Ybar_CC2 = Ybar_CC;
  Ybar_CC3 = Ybar_CC;
  Yvar_CC1 = Yvar_CC;
  Yvar_CC2 = Yvar_CC;
  Yvar_CC3 = Yvar_CC;
}

parameters {
  real<lower = 0, upper = 1>  delta_mu1;       // random mean-discounting parameter specific to EC1
  real<lower = 0, upper = 1>  delta_mu2;       // random mean-discounting parameter specific to EC2
  real<lower = 0, upper = 1>  delta_mu3;       // random mean-discounting parameter specific to EC3
  real  mu_CC1;                                         // population mean of concurrent controls specific to EC1
  real<lower = 0>  sigma_CC1;             // population variance of concurrent controls specific to EC1
  real  mu_CC2;                                         // population mean of concurrent controls specific to EC2
  real<lower = 0>  sigma_CC2;             // population variance of concurrent controls specific to EC2
  real  mu_CC3;                                         // population mean of concurrent controls specific to EC3
  real<lower = 0>  sigma_CC3;             // population variance of concurrent controls specific to EC3
}

transformed parameters{
  real<lower = 0>  sigma_CC_t1;
  real<lower = 0>  sigma_CC_t2;
  real<lower = 0>  sigma_CC_t3;
  
  sigma_CC_t1 = (N_CC1-1)*Yvar_CC1/sigma_CC1;
  sigma_CC_t2 = (N_CC2-1)*Yvar_CC2/sigma_CC2;
  sigma_CC_t3 = (N_CC3-1)*Yvar_CC3/sigma_CC3;
}

model {
  
  // Uniform distribution as the initial prior
  delta_mu1 ~ beta(1, 1);
  delta_mu2 ~ beta(1, 1);
  delta_mu3 ~ beta(1, 1);
  
  // cMCBpP2 constructed from each external control arm with Jeffreys' prior
  target += -1.5*log(sigma_CC1);  
  target += normal_lpdf(mu_CC1|Ybar_EC1, sqrt(sigma_CC1/(N_EC1*delta_mu1)));
  target += -1.5*log(sigma_CC2);  
  target += normal_lpdf(mu_CC2|Ybar_EC2, sqrt(sigma_CC2/(N_EC2*delta_mu2)));
  target += -1.5*log(sigma_CC3);  
  target += normal_lpdf(mu_CC3|Ybar_EC3, sqrt(sigma_CC3/(N_EC3*delta_mu3)));
  
  //likelihood
  target += chi_square_lpdf(sigma_CC_t1|N_CC1-1);
  Ybar_CC1 ~ normal(mu_CC1, sqrt(sigma_CC1/N_CC1));
  target += chi_square_lpdf(sigma_CC_t2|N_CC2-1);
  Ybar_CC2 ~ normal(mu_CC2, sqrt(sigma_CC2/N_CC2));
  target += chi_square_lpdf(sigma_CC_t3|N_CC3-1);
  Ybar_CC3 ~ normal(mu_CC3, sqrt(sigma_CC3/N_CC3));
  
}

generated quantities{

  real<lower = 0>  part1;
  real<lower = 0>  part2;
  real<lower = 0>  part3;
  real<lower = 0>  part4;
  real<lower = 0>  part5;
  real<lower = 0>  N1;
  real<lower = 0>  N2;
  real avr_mu;
  real  mu_CC;                                           // population mean of concurrent controls
  real<lower = 0>  sigma_CC;               // population variance of concurrent controls
  
  N1 = delta_mu1*N_EC1+delta_mu2*N_EC2+delta_mu3*N_EC3+N_CC;
  part1 = Yvar_CC*(N_CC-1);
  part2 = (delta_mu2*N_EC2*delta_mu1*N_EC1)*(Ybar_EC2-Ybar_EC1)^2/(delta_mu1*N_EC1+delta_mu2*N_EC2+delta_mu3*N_EC3);
  part3 = (delta_mu2*N_EC2*delta_mu3*N_EC3)*(Ybar_EC2-Ybar_EC3)^2/(delta_mu1*N_EC1+delta_mu2*N_EC2+delta_mu3*N_EC3);
  part4 = (delta_mu1*N_EC1*delta_mu3*N_EC3)*(Ybar_EC1-Ybar_EC3)^2/(delta_mu1*N_EC1+delta_mu2*N_EC2+delta_mu3*N_EC3);
  part5 = N_CC*(delta_mu1*N_EC1+delta_mu2*N_EC2+delta_mu3*N_EC3)*(Ybar_CC-(delta_mu1*N_EC1*Ybar_EC1+delta_mu2*N_EC2*Ybar_EC2+delta_mu3*N_EC3*Ybar_EC3)/(delta_mu1*N_EC1+delta_mu2*N_EC2+delta_mu3*N_EC3))^2/(N_CC+delta_mu1*N_EC1+delta_mu2*N_EC2+delta_mu3*N_EC3);
  avr_mu = (delta_mu1*N_EC1*Ybar_EC1+delta_mu2*N_EC2*Ybar_EC2+delta_mu3*N_EC3*Ybar_EC3+N_CC*Ybar_CC)/N1;
  N2 = N_CC + 3;
  
  // posterior distribution
  sigma_CC = inv_gamma_rng(N2/2,(part1+part2+part3+part4+part5)/2);
  mu_CC = normal_rng(avr_mu,sqrt(sigma_CC/N1));
  
}


