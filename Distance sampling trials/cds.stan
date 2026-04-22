functions {
  //Define half-normal detection function, to pass to integrate_1d
  // Function arguments here are prescribed by what integrate_1d expects
  // x contains the distance, theta[1] contains the sigma
  // Other arguments not used at present - see help for integrate_1d for more
  real gx_hn(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i){
    return exp(-x ^ 2 / (2 * theta[1] ^ 2));    
  }
}
data {
  int n; //Number of observations
  array[n] real x; //Vector of observed distances
  real w; //Truncation distance
  real L; //Line length
  real<lower = 0.0> hyper_sigma; //standard deviation of half-normal prior on sigma
}
transformed data {
  array[0] real x_r; //Required for integrate_1d - not used
  array[0] int x_i; //Required for integrate_1d - not used
  real er2 = n / (2 * L); //Encounter rate (/2) - used in calculating D
}
parameters {
   real<lower=0.0> sigma; 
}
transformed parameters {
  //Not sure why sigma is in braces - can check out
  real esw_hn = integrate_1d(gx_hn, 0.0, w, { sigma }, x_r, x_i, 1e-8);
}
model {
  //Priors
  sigma ~ normal(0, hyper_sigma); 
  
  //Likelihood
  for(i in 1:n) target +=  (-x[i] ^ 2 / (2 * sigma ^ 2)) - log(esw_hn);
}
generated quantities {
  real p = esw_hn / w;
  //Implement the formula D = n / (2 * L * esw) with n/2L pre-computed as er2 
  real D = er2 / esw_hn; 
}

