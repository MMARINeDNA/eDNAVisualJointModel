// =============================================================================
// visual_functions.stan
//
// Line-transect (visual survey) observation model: half-normal or hazard-rate
// detection with an optional group-size covariate, a group-size distribution, the
// population-average effective strip half-width (ESW), and the size-bias
// correction. Ported (as functions) from the validated
// distance/distance_hn_dens_v4.1.stan, generalised so the detection scale and
// ESW can differ between observations (sea-state covariate): callers pass a
// per-detection base log sigma and per-detection / per-segment ESW. With one
// common value these reduce exactly to v4.1. Included inside `functions { }`:
//
//   functions {
//   #include visual_functions.stan
//   }
//
// Detection key (det_key):
//   0 = half-normal   g(x) = exp(-x^2 / (2 sigma^2))
//   1 = hazard-rate   g(x) = 1 - exp(-(sigma / x)^b)    (ESW by a midpoint rule)
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

  // Hazard-rate ESW  int_0^w [1 - exp(-(sigma/x)^b)] dx  by the midpoint rule on
  // Q cells (log_xq = log of the cell midpoints, built in transformed data).
  // The exponent is capped so exp() cannot overflow; g is 1 there anyway.
  real hr_esw(real log_sigma_c, real b, real w, vector log_xq) {
    vector[rows(log_xq)] t = exp(fmin(b * (fmin(log_sigma_c, 8.0) - log_xq), 50.0));
    return w * mean(-expm1(-t));
  }
  real hr_log_g(real x, real log_sigma_c, real b) {
    real t = exp(fmin(b * (fmin(log_sigma_c, 8.0) - log(fmax(x, 1e-6))), 50.0));
    return log(fmax(-expm1(-t), 1e-300));
  }

  // ESW at a given log scale for the chosen key.
  real det_esw(real log_sigma_c, int det_key, real b, real w, vector log_xq) {
    if (det_key == 0) return hn_esw(hn_sigma(log_sigma_c), w);
    return hr_esw(log_sigma_c, b, w, log_xq);
  }

  // Population-average ESW: E_pop[ESW(sigma(s))] over the group-size pmf.
  // Without a size covariate sigma is constant and this is just det_esw().
  //
  // The sum over group sizes is done on kbin bins (each size k <= 30 its own
  // bin, larger sizes in log bins of width 10%): the bin's pmf mass at its
  // pmf-weighted mean size. ESW(s) is smooth in s, so this matches the full
  // sum to ~1e-5 while cutting the ESW evaluations from S_max to ~70.
  real esw_population(vector pk, int use_size_covar, int det_key, real log_sigma,
                      real beta_size, real s_centre, real b, real w,
                      vector log_xq, array[] int kbin, int nkb) {
    if (use_size_covar == 0) {
      return det_esw(log_sigma, det_key, b, w, log_xq);
    }
    int K = rows(pk);
    vector[nkb] pb = rep_vector(0, nkb);
    vector[nkb] kb = rep_vector(0, nkb);
    for (k in 1:K) {
      pb[kbin[k]] += pk[k];
      kb[kbin[k]] += pk[k] * k;
    }
    vector[nkb] dk = kb ./ fmax(pb, 1e-300) - s_centre;   // mean size in bin - s_centre
    if (det_key == 0) {
      vector[nkb] sg = exp(fmin(log_sigma + beta_size * dk, 30.0));
      return dot_product(pb, sg * sqrt(pi() / 2) .* erf(w ./ (sqrt(2) * sg)));
    }
    real esw = 0;
    for (j in 1:nkb) {
      esw += pb[j] * det_esw(log_sigma + beta_size * dk[j], det_key, b, w, log_xq);
    }
    return esw;
  }

  // Per-detection log-likelihood: perpendicular distance + group size with
  // the size-bias correction. p(x, s | detected) = g(x | s) f_pop(s) / esw_pop,
  // i.e. log g(x_i) + log f_pop(s_i) - log esw_pop_i (the per-detection ESW
  // terms of v4.1's derivation cancel).
  // log_sigma[i]: detection i's log sigma at s = s_centre (includes any
  // sea-state term); esw_pop[i]: the population-average ESW at that log sigma.
  vector lt_detection_loglik(vector x, array[] int s, real w,
                             int use_size_covar, int det_key, vector log_sigma,
                             real beta_size, real s_centre, real b,
                             int group_size_dist, real mu_s, real phi_s,
                             real mu_log_s, real sigma_log_s, vector esw_pop) {
    int n = rows(x);
    vector[n] ll;
    real log_p0_trunc = group_size_dist == 0
                        ? log1m(neg_binomial_2_cdf(0 | mu_s, phi_s)) : 0;
    for (i in 1:n) {
      real log_sigma_i = use_size_covar == 1
                         ? log_sigma[i] + beta_size * (s[i] - s_centre)
                         : log_sigma[i];
      real log_fs  = group_size_dist == 0
                     ? neg_binomial_2_lpmf(s[i] | mu_s, phi_s) - log_p0_trunc
                     : lognormal_lpdf(s[i] | mu_log_s, sigma_log_s);
      real log_g   = det_key == 0
                     ? -square(x[i]) / (2.0 * square(hn_sigma(log_sigma_i)))
                     : hr_log_g(x[i], log_sigma_i, b);
      ll[i] = log_g + log_fs - log(esw_pop[i]);
    }
    return ll;
  }

  // Per-segment Poisson encounter log-likelihood. log_lambda_animals is the
  // log ANIMAL density (animals / km^2, the field shared with eDNA);
  // expected detected groups = lambda_animals / E[s] * 2 * L * esw_pop[j]
  // (esw_pop[j]: population-average ESW under segment j's conditions).
  vector lt_segment_loglik(array[] int seg_count, vector seg_l,
                           vector log_lambda_animals, real log_mean_group_size,
                           vector esw_pop) {
    int J = size(seg_count);
    vector[J] ll;
    for (j in 1:J) {
      ll[j] = poisson_log_lpmf(seg_count[j] | log_lambda_animals[j] - log_mean_group_size
                                              + log(2 * seg_l[j] * esw_pop[j]));
    }
    return ll;
  }

  // Negative-binomial version (clustered sightings): the same mean, variance
  // mu + mu^2 / phi_nb. A gamma(phi_nb, phi_nb) multiplier on the groups in a
  // segment's strip, thinned by detection, gives exactly this.
  vector lt_segment_loglik_nb(array[] int seg_count, vector seg_l,
                              vector log_lambda_animals, real log_mean_group_size,
                              vector esw_pop, real phi_nb) {
    int J = size(seg_count);
    vector[J] ll;
    for (j in 1:J) {
      ll[j] = neg_binomial_2_log_lpmf(seg_count[j] | log_lambda_animals[j] - log_mean_group_size
                                                     + log(2 * seg_l[j] * esw_pop[j]), phi_nb);
    }
    return ll;
  }
