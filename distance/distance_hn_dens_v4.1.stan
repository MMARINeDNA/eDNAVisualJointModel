// =============================================================================
// distance_hn_dens_v4.1.stan
//
// Single-species, non-spatial line-transect distance-sampling model with:
//   * Half-normal detection function with optional group-size covariate
//     on log(sigma):  log(sigma_i) = log_sigma + beta_size * (s_i - s_centre)
//   * Zero-truncated negative-binomial group-size distribution
//   * Poisson encounter rate per segment
//
// Generalisation of distance/distance_hn_dens.stan. Reduces to that
// model when use_size_covar = 0 (and given matching priors).
//
// New in v4.1:
//   * log_lambda_s prior parameters now passed as data (was hardcoded
//     normal(-8, 2) in the original).
//   * mu_s and phi_s prior parameters passed as data (was hardcoded
//     gamma(2, 0.02) and gamma(1, 0.1)). Allows species-specific
//     priors - critical for humpbacks where mu_s ~ 2 vs PWSD ~ 50.
//   * Optional group-size covariate on the detection function:
//        sigma_i = exp(log_sigma + beta_size * (s_i - s_centre))
//     with a Normal prior on beta_size also passed as data. Set
//     use_size_covar = 0 to disable (beta_size still sampled, but
//     tightly pinned by its prior so it doesn't matter).
//   * For the encounter-rate Poisson, the representative sigma is
//     evaluated at the *population* mean group size (= mu_s, the NB
//     parameter), so size-biased detection is accounted for.
// =============================================================================

data {
  int<lower=1> n;                          // total detected groups
  vector<lower=0>[n] x;                    // observed perpendicular distances
  array[n] int<lower=1> s;                 // observed group sizes
  real<lower=0> w;                         // truncation distance

  // Detection-function priors
  real log_sigma_prior_mean;
  real<lower=0> log_sigma_prior_sd;

  // Group-size covariate
  int<lower=0, upper=1> use_size_covar;
  real s_centre;                           // centring constant for s
  real beta_size_prior_mean;
  real<lower=0> beta_size_prior_sd;

  // Group-size NB priors (now species-specific)
  real<lower=0> mu_s_prior_shape;
  real<lower=0> mu_s_prior_rate;
  real<lower=0> phi_s_prior_shape;
  real<lower=0> phi_s_prior_rate;

  // Group-density prior (now species-specific)
  real log_lambda_s_prior_mean;
  real<lower=0> log_lambda_s_prior_sd;

  // Encounter-rate data
  int<lower=1> n_seg;
  array[n_seg] int<lower=0> seg_count;
  vector<lower=0>[n_seg] seg_l;
}

transformed data {
  vector[n] s_c = to_vector(s) - s_centre;
}

parameters {
  real log_sigma;             // log(sigma) at s = s_centre
  real beta_size;             // log_sigma slope on (s - s_centre)
  real<lower=0> mu_s;
  real<lower=0> phi_s;
  real log_lambda_s;
}

transformed parameters {
  // Per-observation HN scale and ESW.
  vector<lower=0>[n] sigma_i;
  vector<lower=0>[n] esw_i;

  for (i in 1:n) {
    if (use_size_covar == 1)
      sigma_i[i] = exp(log_sigma + beta_size * s_c[i]);
    else
      sigma_i[i] = exp(log_sigma);
    esw_i[i] = sigma_i[i] * sqrt(pi() / 2)
               * erf(w / (sqrt(2) * sigma_i[i]));
  }

  // Representative sigma / ESW for the encounter-rate Poisson, at the
  // POPULATION mean group size (= mu_s, the NB parameter). Approximate
  // way to fold the size distribution into the encounter rate while
  // still using a single ESW per segment.
  real<lower=0> sigma_rep;
  if (use_size_covar == 1)
    sigma_rep = exp(log_sigma + beta_size * (mu_s - s_centre));
  else
    sigma_rep = exp(log_sigma);
  real<lower=0> esw_rep = sigma_rep * sqrt(pi() / 2)
                          * erf(w / (sqrt(2) * sigma_rep));

  real<lower=0> lambda_s = exp(log_lambda_s);
}

model {
  // Priors
  log_sigma    ~ normal(log_sigma_prior_mean,    log_sigma_prior_sd);
  beta_size    ~ normal(beta_size_prior_mean,    beta_size_prior_sd);
  mu_s         ~ gamma(mu_s_prior_shape,         mu_s_prior_rate);
  phi_s        ~ gamma(phi_s_prior_shape,        phi_s_prior_rate);
  log_lambda_s ~ normal(log_lambda_s_prior_mean, log_lambda_s_prior_sd);

  // HN detection-function likelihood (per-observation)
  for (i in 1:n) {
    target += -square(x[i]) / (2.0 * square(sigma_i[i])) - log(esw_i[i]);
  }

  // Zero-truncated NB on group sizes
  target += neg_binomial_2_lpmf(s | mu_s, phi_s)
            - n * log1m(neg_binomial_2_cdf(0 | mu_s, phi_s));

  // Encounter rate Poisson per segment, using representative ESW
  target += poisson_lpmf(seg_count | lambda_s * seg_l * 2 * esw_rep);
}

generated quantities {
  // Average detection probability in the strip at population mean size
  real<lower=0, upper=1> p = esw_rep / w;
  // Group density (groups/km^2) and animal density (animals/km^2).
  // Note the original distance_hn_dens.stan reported D_s in groups per
  // 1000 area units; we report in animals/km^2 directly.
  real<lower=0> D_animals = lambda_s * mu_s;
}
