//
//  Modified Power Prior
//  To leverage external Gaussian endpoint with unknown variance (MPP1) from three external control arms
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

parameters {
  real<lower = 0, upper = 1>  delta1;  // random power parameter specific to EC1
  real<lower = 0, upper = 1>  delta2;  // random power parameter specific to EC2
  real<lower = 0, upper = 1>  delta3;  // random power parameter specific to EC3
  real  mu_CC;                                            // population mean of concurrent controls
  real<lower = 0>  sigma_CC;                // population variance of concurrent controls
}

transformed parameters{

  real<lower = 0>  part1;
  real<lower = 0>  part2;
  real<lower = 0>  part3;
  real<lower = 0>  part4;
  real<lower = 0>  N;
  real avr_mu;
  real<lower = 0>  sigma_CC_t;
  
  N          = delta1*N_EC1 + delta2*N_EC2 + delta3*N_EC3;
  part1      = delta1*Yvar_EC1*(N_EC1-1) + delta2*Yvar_EC2*(N_EC2-1) + delta3*Yvar_EC3*(N_EC3-1);
  part2      = (delta2*N_EC2*delta1*N_EC1)*(Ybar_EC2-Ybar_EC1)^2/N;
  part3      = (delta2*N_EC2*delta3*N_EC3)*(Ybar_EC2-Ybar_EC3)^2/N;
  part4      = (delta1*N_EC1*delta3*N_EC3)*(Ybar_EC1-Ybar_EC3)^2/N;
  avr_mu     = (delta1*N_EC1*Ybar_EC1+delta2*N_EC2*Ybar_EC2+delta3*N_EC3*Ybar_EC3)/N;
  sigma_CC_t = (N_CC-1)*Yvar_CC/sigma_CC;
}


model {
  
  // Uniform distribution as the initial prior
  delta1 ~ beta(1, 1);
  delta2 ~ beta(1, 1);
  delta3 ~ beta(1, 1);
  
  // MPP1 constructed using all external controls and the Jeffreys' prior
  target += inv_gamma_lpdf(sigma_CC|N/2, (part1+part2+part3+part4)/2);
  target += normal_lpdf(mu_CC|avr_mu, sqrt(sigma_CC/N));
  
  //likelihood
  target += chi_square_lpdf(sigma_CC_t|N_CC-1);
  Ybar_CC ~ normal(mu_CC, sqrt(sigma_CC/N_CC));
  
}

