data {
  int<lower=1> n; //Number of observations
  vector<lower=0>[n] x; //Vector of observed distances
  array[n] int<lower=1> s; //Vector of observed group sizes
  real<lower=0> w; //Truncation distance
  real log_sigma_prior_mean; //Mean of normal prior on log sigma
  real<lower=0> log_sigma_prior_sd; //SD of normal prior on log sigma
  int<lower=1> n_seg; //Number of segments
  array[n_seg] int<lower=0> seg_count; //Number of detections per segment 
  vector<lower=0>[n_seg] seg_l; //Line length per segment
}

parameters {
  real log_sigma; //Log of half normal sigma parameter
  real<lower=0> mu_s; //Population mean group size
  real<lower=0> phi_s; //Overdispersion in neg bin group size distribution
  real log_lambda; //Log of group density
}

transformed parameters {
  real<lower=0> sigma = exp(log_sigma); //Define sigma
  //Calculate effective strip width - can be done analytically in this case
  // as the detection function is a half normal with no adjustments
  real<lower=0> esw_hn = sigma * sqrt(pi() / 2) * erf(w / (sqrt(2) * sigma));
  real<lower=0> lambda = exp(log_lambda); //Define lambda: group density
}

model {
  //Priors
  log_sigma ~ normal(log_sigma_prior_mean, log_sigma_prior_sd); 
  mu_s ~ gamma(2, 0.02); //Weakly informative prior, apparently
  phi_s ~ gamma(1, 0.1); //Allows strong overdispersion
  //Value close to average seg_count/(2lw) - need to revisit to make general
  log_lambda ~ normal(-8, 2); 

  //Likelihood - vectorized
  //HN detection function
  target += -sum(square(x)) / (2.0 * square(sigma)) - n * log(esw_hn);
  //Zero truncated negative binomial on group size  -- implemented
  // using neg bin lpmf and subtracting off probability of zero
  target += neg_binomial_2_lpmf(s | mu_s, phi_s)
            - n * log1m(neg_binomial_2_cdf(0 | mu_s, phi_s));
  //Poisson likelihood for encounter rate
  target += poisson_lpmf(seg_count | lambda * seg_l * 2 * esw_hn);
}

generated quantities {
  //Average prob of detection between 0 and w
  real<lower=0, upper=1> p = esw_hn / w;
  //Density of groups in ha
  real<lower=0> D_s = lambda * 100;
  //Density of animals in ha
  real<lower=0> D = D_s * mu_s;
  

//  // Posterior predictive draws for group size
//  array[n] int s_rep;
//  for (i in 1:n) {
//    int draw;
//    draw = neg_binomial_2_rng(mu_s, phi_s);
//    // Enforce zero truncation
//    while (draw == 0)
//      draw = neg_binomial_2_rng(mu_s, phi_s);
//    s_rep[i] = draw;
//  }
}
