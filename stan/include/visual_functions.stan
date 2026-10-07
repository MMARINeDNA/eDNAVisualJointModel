// =============================================================================
// visual_functions.stan
//
// Line-transect (visual survey) observation model: half-normal detection with
// an optional group-size covariate, a group-size distribution, the
// population-average effective strip half-width (ESW), and the size-bias
// correction. Ported verbatim (as functions) from the validated
// distance/distance_hn_dens_v4.1.stan. Included inside `functions { }`:
//
//   functions {
//   #include visual_functions.stan
//   }
//
// Group-size distribution (group_size_dist):
//   0 = zero-truncated NB(mu_s, phi_s)
//   1 = log-normal(mu_log_s, sigma_log_s), rounded to integers k = 1..S_max
// =============================================================================

  // Half-normal ESW for scale sigma, truncation w:  int_0^w exp(-x^2 / 2 sigma^2) dx.
  // log(sigma) is clamped at 30 (sigma ~ 1e13 km): far beyond any truncation
  // distance, where ESW is flat at w, so the clamp never binds in a plausible
  // posterior region. It only stops exp() overflowing during warmup when
  // beta_size wanders (log_sigma + beta_size * (k - s_centre) for k up to S_max).
  real hn_sigma(real log_sigma) {
    return exp(fmin(log_sigma, 30.0));
  }
  real hn_esw(real sigma, real w) {
    return sigma * sqrt(pi() / 2) * erf(w / (sqrt(2) * sigma));
  }

  // Population group-size pmf on k = 1..S_max, normalised to sum to 1.
  vector group_size_pmf(int S_max, int group_size_dist, real mu_s, real phi_s,
                        real mu_log_s, real sigma_log_s) {
    vector[S_max] pk;
    for (k in 1:S_max) {
      if (group_size_dist == 0) {
        pk[k] = exp(neg_binomial_2_lpmf(k | mu_s, phi_s));
      } else {
        real lo = log(k - 0.5);   // log(0.5) for k = 1
        real hi = log(k + 0.5);
        pk[k] = Phi((hi - mu_log_s) / sigma_log_s)
                - Phi((lo - mu_log_s) / sigma_log_s);
      }
    }
    return pk / sum(pk);
  }

  // Population-average ESW: E_pop[ESW(sigma(s))] over the group-size pmf.
  // Without a size covariate sigma is constant and this is just hn_esw().
  real esw_population(vector pk, int use_size_covar, real log_sigma,
                      real beta_size, real s_centre, real w) {
    if (use_size_covar == 0) {
      return hn_esw(hn_sigma(log_sigma), w);
    }
    real esw = 0;
    for (k in 1:rows(pk)) {
      esw += pk[k] * hn_esw(hn_sigma(log_sigma + beta_size * (k - s_centre)), w);
    }
    return esw;
  }

  // Per-detection log-likelihood: perpendicular distance + group size with
  // the size-bias correction.
  //   distance:   -x^2 / (2 sigma_i^2) - log(esw_i)           (truncated HN pdf, up to -log w)
  //   group size: log f_pop(s) + log(esw_i) - log(esw_pop)     (f(s | detected))
  // The log(esw_i) terms cancel; both are kept to mirror v4.1's derivation.
  vector lt_detection_loglik(vector x, array[] int s, real w,
                             int use_size_covar, real log_sigma,
                             real beta_size, real s_centre,
                             int group_size_dist, real mu_s, real phi_s,
                             real mu_log_s, real sigma_log_s, real esw_pop) {
    int n = rows(x);
    vector[n] ll;
    real log_p0_trunc = group_size_dist == 0
                        ? log1m(neg_binomial_2_cdf(0 | mu_s, phi_s)) : 0;
    for (i in 1:n) {
      real sigma_i = hn_sigma(use_size_covar == 1
                              ? log_sigma + beta_size * (s[i] - s_centre)
                              : log_sigma);
      real esw_i   = hn_esw(sigma_i, w);
      real log_fs  = group_size_dist == 0
                     ? neg_binomial_2_lpmf(s[i] | mu_s, phi_s) - log_p0_trunc
                     : lognormal_lpdf(s[i] | mu_log_s, sigma_log_s);
      ll[i] = -square(x[i]) / (2.0 * square(sigma_i)) - log(esw_i)
              + log_fs + log(esw_i) - log(esw_pop);
    }
    return ll;
  }

  // Per-segment Poisson encounter log-likelihood. log_lambda_animals is the
  // log ANIMAL density (animals / km^2, the field shared with eDNA);
  // expected detected groups = lambda_animals / E[s] * 2 * L * esw_pop.
  vector lt_segment_loglik(array[] int seg_count, vector seg_l,
                           vector log_lambda_animals, real log_mean_group_size,
                           real esw_pop) {
    int J = size(seg_count);
    vector[J] ll;
    for (j in 1:J) {
      ll[j] = poisson_log_lpmf(seg_count[j] | log_lambda_animals[j] - log_mean_group_size
                                              + log(2 * seg_l[j] * esw_pop));
    }
    return ll;
  }
