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

  // Upper integration limit for the population-average ESW (zero-truncated NB
  // population-mean ESW). Should be large enough to cover the right tail of
  // NB(mu_s, phi_s) under any plausible posterior values; e.g. for PWSD
  // (mu_s ~ 50, phi_s small) S_max = 500-1000 is safe.
  int<lower=1> S_max;

  // Group-size distribution choice:
  //   0 = zero-truncated NB(mu_s, phi_s)  (Option default; works for
  //        moderately overdispersed counts, e.g. humpback)
  //   1 = log-normal on s (Option A; better for heavy right tails like
  //        PWSD where the empirical mean ~ 43, max ~ 450)
  int<lower=0, upper=1> group_size_dist;

  // Group-size NB priors (used when group_size_dist == 0)
  real<lower=0> mu_s_prior_shape;
  real<lower=0> mu_s_prior_rate;
  real<lower=0> phi_s_prior_shape;
  real<lower=0> phi_s_prior_rate;

  // Group-size log-normal priors (used when group_size_dist == 1)
  real mu_log_prior_mean;
  real<lower=0> mu_log_prior_sd;
  real<lower=0> sigma_log_prior_shape;
  real<lower=0> sigma_log_prior_rate;

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
  // NB group-size params (used iff group_size_dist == 0; otherwise
  // sampled but unused, just sitting at their prior)
  real<lower=0> mu_s;
  real<lower=0> phi_s;
  // Log-normal group-size params (used iff group_size_dist == 1)
  real           mu_log_s;
  real<lower=0>  sigma_log_s;
  real log_lambda_s;
}

transformed parameters {
  // Per-observation HN scale and ESW. The `log(sigma)` argument to `exp`
  // is clamped at 30 (≈ 1e13, far beyond any plausible `w`) so the chain
  // can't blow `sigma` to Inf during warmup wandering. With `beta_size`
  // proposed at, say, 0.7 and PWSD's `s_centre = 43, S_max = 1000`,
  // `log_sigma + 0.7 * (k − 43)` reaches ~670 inside the integration
  // loop below, which would overflow `exp` to Inf; then `Inf * erf(small)`
  // is NaN and the `<lower=0>` check on `esw_pop` rejects every step.
  // `esw = sigma * sqrt(π/2) * erf(w/(sigma√2))` converges cleanly to
  // `w` as sigma → ∞, so the clamp doesn't change the answer in any
  // posterior region — it only prevents the numerical pathology.
  vector<lower=0>[n] sigma_i;
  vector<lower=0>[n] esw_i;

  for (i in 1:n) {
    real log_sigma_i;
    if (use_size_covar == 1)
      log_sigma_i = log_sigma + beta_size * s_c[i];
    else
      log_sigma_i = log_sigma;
    sigma_i[i] = exp(fmin(log_sigma_i, 30.0));
    esw_i[i] = sigma_i[i] * sqrt(pi() / 2)
               * erf(w / (sqrt(2) * sigma_i[i]));
  }

  // Representative sigma / ESW at the POPULATION MEAN group size
  // (= mu_s, the NB parameter). Kept for reporting; not used in the
  // encounter-rate likelihood.
  real<lower=0> sigma_rep;
  {
    real log_sigma_rep;
    if (use_size_covar == 1)
      log_sigma_rep = log_sigma + beta_size * (mu_s - s_centre);
    else
      log_sigma_rep = log_sigma;
    sigma_rep = exp(fmin(log_sigma_rep, 30.0));
  }
  real<lower=0> esw_rep = sigma_rep * sqrt(pi() / 2)
                          * erf(w / (sqrt(2) * sigma_rep));

  // Population-average ESW: numerical integration of
  // ESW(sigma(s)) over the population group-size distribution. With
  // sigma(s) = exp(log_sigma + beta * (s - s_centre)) convex in s,
  // E[ESW(sigma(s))] > ESW(sigma(E[s])) by Jensen's inequality, so using esw_rep in
  // the encounter-rate Poisson under-estimates the true mean ESW and
  // inflates lambda_s. esw_pop fixes that. Sum runs over k = 1, ...,
  // S_max with pmf P(k) under the chosen group-size distribution.
  // For log-normal, we use the rounded-bin pmf
  //   P(k) = Phi((log(k+0.5) - mu_log)/sigma_log)
  //          - Phi((log(k-0.5) - mu_log)/sigma_log).
  real<lower=0> esw_pop;
  if (use_size_covar == 1) {
    real esw_sum = 0;
    real p_sum   = 0;
    for (k in 1:S_max) {
      real pk;
      if (group_size_dist == 0) {
        pk = exp(neg_binomial_2_lpmf(k | mu_s, phi_s));
      } else {
        real lo = log(k - 0.5);   // log(0.5) for k = 1
        real hi = log(k + 0.5);
        pk = Phi((hi - mu_log_s) / sigma_log_s)
             - Phi((lo - mu_log_s) / sigma_log_s);
      }
      real sigma_k = exp(fmin(log_sigma + beta_size * (k - s_centre), 30.0));
      real esw_k   = sigma_k * sqrt(pi() / 2)
                     * erf(w / (sqrt(2) * sigma_k));
      esw_sum += pk * esw_k;
      p_sum   += pk;
    }
    esw_pop = esw_sum / p_sum;
  } else {
    esw_pop = esw_rep;
  }

  real<lower=0> lambda_s = exp(log_lambda_s);
}

model {
  // Priors. Both NB and log-normal group-size priors are always
  // sampled; only the one matching `group_size_dist` enters the data
  // likelihood. The other gets its prior but doesn't affect inference.
  log_sigma    ~ normal(log_sigma_prior_mean,    log_sigma_prior_sd);
  beta_size    ~ normal(beta_size_prior_mean,    beta_size_prior_sd);
  mu_s         ~ gamma(mu_s_prior_shape,         mu_s_prior_rate);
  phi_s        ~ gamma(phi_s_prior_shape,        phi_s_prior_rate);
  mu_log_s     ~ normal(mu_log_prior_mean,       mu_log_prior_sd);
  sigma_log_s  ~ gamma(sigma_log_prior_shape,    sigma_log_prior_rate);
  log_lambda_s ~ normal(log_lambda_s_prior_mean, log_lambda_s_prior_sd);

  // HN detection-function likelihood (per-observation)
  for (i in 1:n) {
    target += -square(x[i]) / (2.0 * square(sigma_i[i])) - log(esw_i[i]);
  }

  // Group-size likelihood, plus the size-biased detection correction.
  // Without the correction term, observed sizes are treated as iid
  // draws from the population, which is wrong: when detection
  // probability depends on size (use_size_covar = 1), larger groups
  // are over-represented in the detected sample. The correction
  // conditions on detection:
  //   f(s | detected) propto f_pop(s) * p_det(s)
  //   log f(s|det) = log f_pop(s) + log(p_det(s)) - log(E_pop[p_det])
  //                = log f_pop(s) + log(esw_i) - log(esw_pop)
  // (The -log(w) constants cancel.) When use_size_covar = 0, esw_i is
  // constant and the correction term is identically 0.
  if (group_size_dist == 0) {
    // Zero-truncated NB
    target += neg_binomial_2_lpmf(s | mu_s, phi_s)
              - n * log1m(neg_binomial_2_cdf(0 | mu_s, phi_s));
  } else {
    // Log-normal on s (continuous-approximation Option A). Stan's
    // lognormal_lpdf gives the density on the s scale directly.
    target += lognormal_lpdf(s | mu_log_s, sigma_log_s);
  }
  target += sum(log(esw_i)) - n * log(esw_pop);

  // Encounter rate Poisson per segment, using the population-average
  // ESW (averaged over zero-truncated NB).
  target += poisson_lpmf(seg_count | lambda_s * seg_l * 2 * esw_pop);
}

generated quantities {
  // Mean detection probability in the strip, averaged over the
  // population group-size distribution (population-average).
  real<lower=0, upper=1> p = esw_pop / w;
  // Same quantity at the population mean size (for comparison).
  real<lower=0, upper=1> p_at_mu = esw_rep / w;
  // Group density (groups/km^2) and animal density (animals/km^2).
  // Note the original distance_hn_dens.stan reported D_s in groups per
  // 1000 area units; we report in animals/km^2 directly.
  real<lower=0> D_animals = lambda_s * mu_s;
}
