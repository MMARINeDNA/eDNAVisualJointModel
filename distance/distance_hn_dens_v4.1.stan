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

  // Upper integration limit for the Jensen-corrected ESW (zero-truncated NB
  // population-mean ESW). Should be large enough to cover the right tail of
  // NB(mu_s, phi_s) under any plausible posterior values; e.g. for PWSD
  // (mu_s ~ 50, phi_s small) S_max = 500-1000 is safe.
  int<lower=1> S_max;

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

  // Representative sigma / ESW at the POPULATION MEAN group size
  // (= mu_s, the NB parameter). Kept for reporting; not used in the
  // encounter-rate likelihood.
  real<lower=0> sigma_rep;
  if (use_size_covar == 1)
    sigma_rep = exp(log_sigma + beta_size * (mu_s - s_centre));
  else
    sigma_rep = exp(log_sigma);
  real<lower=0> esw_rep = sigma_rep * sqrt(pi() / 2)
                          * erf(w / (sqrt(2) * sigma_rep));

  // Jensen-corrected population ESW: numerical integration of
  // ESW(sigma(s)) over the zero-truncated NB(mu_s, phi_s)
  // distribution. With sigma(s) = exp(log_sigma + beta * (s - s_centre))
  // convex in s, E[ESW(sigma(s))] > ESW(sigma(mu_s)) by Jensen, so
  // using esw_rep above in the encounter rate Poisson under-estimates
  // the true mean ESW and inflates the lambda_s posterior to compensate.
  // esw_pop fixes that.
  real<lower=0> esw_pop;
  if (use_size_covar == 1) {
    real esw_sum = 0;
    real p_sum   = 0;
    for (k in 1:S_max) {
      real pk      = exp(neg_binomial_2_lpmf(k | mu_s, phi_s));
      real sigma_k = exp(log_sigma + beta_size * (k - s_centre));
      real esw_k   = sigma_k * sqrt(pi() / 2)
                     * erf(w / (sqrt(2) * sigma_k));
      esw_sum += pk * esw_k;
      p_sum   += pk;
    }
    // Renormalise (pmf truncated at S_max; also drops the s = 0 mass
    // implicitly since we start at k = 1).
    esw_pop = esw_sum / p_sum;
  } else {
    // Without the covariate, sigma is constant in s -> integration
    // collapses to ESW at sigma_0.
    esw_pop = esw_rep;
  }

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

  // Zero-truncated NB on group sizes, with size-biased detection
  // correction. Without the correction term, observed sizes are
  // treated as iid draws from the population NB, which is wrong:
  // when detection probability depends on size (use_size_covar = 1),
  // larger groups are over-represented in the detected sample, and
  // fitting NB to the observed sample biases mu_s upward. The
  // correction conditions on detection:
  //   f(s | detected) propto f_pop(s) * p_det(s)
  //   log f(s|det) = log NB(s|mu, phi) - log(1 - P(s=0))
  //                 + log(p_det(s)) - log(E_pop[p_det])
  //                = ... + log(esw_i) - log(esw_pop)
  // The constant -log(w) cancels. When use_size_covar = 0, esw_i is
  // the same for all i and the correction is 0.
  target += neg_binomial_2_lpmf(s | mu_s, phi_s)
            - n * log1m(neg_binomial_2_cdf(0 | mu_s, phi_s));
  target += sum(log(esw_i)) - n * log(esw_pop);

  // Encounter rate Poisson per segment, using Jensen-corrected
  // population-mean ESW (averaged over zero-truncated NB).
  target += poisson_lpmf(seg_count | lambda_s * seg_l * 2 * esw_pop);
}

generated quantities {
  // Mean detection probability in the strip, averaged over the
  // population group-size distribution (Jensen-corrected).
  real<lower=0, upper=1> p = esw_pop / w;
  // Same quantity at the population mean size (for comparison).
  real<lower=0, upper=1> p_at_mu = esw_rep / w;
  // Group density (groups/km^2) and animal density (animals/km^2).
  // Note the original distance_hn_dens.stan reported D_s in groups per
  // 1000 area units; we report in animals/km^2 directly.
  real<lower=0> D_animals = lambda_s * mu_s;
}
