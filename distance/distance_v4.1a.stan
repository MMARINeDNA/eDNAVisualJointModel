// =============================================================================
// distance_v4.1a.stan
//
// SPATIAL extension of distance/distance_hn_dens_v4.1.stan. The
// non-spatial model fit a single scalar group density log_lambda_s; this
// model replaces that scalar with a 3-D anisotropic HSGP latent field
// over segment locations, mirroring the spatial structure of
// stan/whale_edna_hsgp_v3.2.stan but with the line-transect observation
// model from distance_hn_dens_v4.1.stan.
//
// Latent field (single species per fit, S = 1):
//   log lambda_groups(x_j, y_j, z_j) = mu_sp + f(x_j, y_j, z_j)
//   f ~ GP(0, K(lx, ly, lz))   approximated via HSGP basis
//
// Observation model (line-transect, cetaceans):
//   * seg_count[j] ~ Poisson(2 * seg_l[j] * esw_pop * lambda_groups[j])
//   * x_i           ~ Half-normal(sigma_i)  truncated at w
//   * group sizes   ~ ZTNB(mu_s, phi_s)  OR  log-normal(mu_log, sigma_log)
//
// All detection-function machinery (population-corrected ESW with optional
// group-size covariate, size-bias correction in the group-size
// likelihood) is identical to distance_hn_dens_v4.1.stan -- the only
// change is that lambda_groups is now per-segment instead of a single
// scalar.
//
// Generated quantities adds posterior-predictive checks for:
//   * lambda_groups / lambda_animals at each segment (the spatial field)
//   * f at each segment (the spatial component, mu_sp removed)
//   * simulated segment counts (PPC for the encounter-rate model)
//   * per-segment + per-detection log-likelihoods (for loo / PSIS)
// =============================================================================

functions {

  // -------------------------------------------------------------------------
  // 1-D Laplacian eigenfunction on [-L, L]   (verbatim from v3.2)
  // -------------------------------------------------------------------------
  real phi1d(real L, int m, real x) {
    return (1.0 / sqrt(L)) * sin(m * pi() * (x + L) / (2.0 * L));
  }

  // -------------------------------------------------------------------------
  // Spectral density of 1-D squared-exponential kernel
  // -------------------------------------------------------------------------
  real spd_se_1d(real alpha, real rho, real w_) {
    return square(alpha) * sqrt(2.0 * pi()) * rho
             * exp(-0.5 * square(rho * w_));
  }

  // -------------------------------------------------------------------------
  // N x M HSGP design matrix (3-D product of 1-D eigenfunctions)
  // -------------------------------------------------------------------------
  matrix hsgp_phi(matrix coords, vector L, array[] int m) {
    int N = rows(coords);
    int M = m[1] * m[2] * m[3];
    matrix[N, M] PHI;
    int col = 1;
    for (j1 in 1:m[1]) {
      for (j2 in 1:m[2]) {
        for (j3 in 1:m[3]) {
          for (n in 1:N) {
            PHI[n, col] =
              phi1d(L[1], j1, coords[n, 1]) *
              phi1d(L[2], j2, coords[n, 2]) *
              phi1d(L[3], j3, coords[n, 3]);
          }
          col += 1;
        }
      }
    }
    return PHI;
  }

  // -------------------------------------------------------------------------
  // Spectral weights (length M) for anisotropic 3-D SE kernel
  // -------------------------------------------------------------------------
  vector hsgp_weights(real sigma_gp, vector l, vector L,
                      vector coord_scale, array[] int m) {
    int M = m[1] * m[2] * m[3];
    vector[M] wt;
    real alpha_dim = pow(sigma_gp, 1.0 / 3.0);
    int col = 1;
    for (j1 in 1:m[1]) {
      for (j2 in 1:m[2]) {
        for (j3 in 1:m[3]) {
          real freq1 = (j1 * pi() / (2.0 * L[1])) / coord_scale[1];
          real freq2 = (j2 * pi() / (2.0 * L[2])) / coord_scale[2];
          real freq3 = (j3 * pi() / (2.0 * L[3])) / coord_scale[3];
          wt[col] = sqrt(
            spd_se_1d(alpha_dim, l[1], freq1) *
            spd_se_1d(alpha_dim, l[2], freq2) *
            spd_se_1d(alpha_dim, l[3], freq3)
          );
          col += 1;
        }
      }
    }
    return wt;
  }

}

// =============================================================================
data {

  // ---- Detections (one row per detected group) ------------------------------
  int<lower=1> n;                          // total detected groups
  vector<lower=0>[n] x;                    // observed perpendicular distances
  array[n]    int<lower=1> s;              // observed group sizes
  array[n]    int<lower=1> det_seg;        // segment index for each detection
  real<lower=0> w;                         // truncation distance

  // ---- Detection-function priors ------------------------------------------
  real log_sigma_prior_mean;
  real<lower=0> log_sigma_prior_sd;
  int<lower=0, upper=1> use_size_covar;
  real s_centre;                           // centring constant for s
  real beta_size_prior_mean;
  real<lower=0> beta_size_prior_sd;

  // ---- Population-corrected ESW integration cap ----------------------------------
  int<lower=1> S_max;

  // ---- Group-size sub-model -----------------------------------------------
  // 0 = zero-truncated NB(mu_s, phi_s);  1 = log-normal(mu_log, sigma_log)
  int<lower=0, upper=1> group_size_dist;

  // NB priors (used iff group_size_dist == 0)
  real<lower=0> mu_s_prior_shape;
  real<lower=0> mu_s_prior_rate;
  real<lower=0> phi_s_prior_shape;
  real<lower=0> phi_s_prior_rate;

  // Log-normal priors (used iff group_size_dist == 1)
  real           mu_log_prior_mean;
  real<lower=0>  mu_log_prior_sd;
  real<lower=0>  sigma_log_prior_shape;
  real<lower=0>  sigma_log_prior_rate;

  // ---- Encounter-rate data (per segment) ----------------------------------
  int<lower=1>                  n_seg;
  array[n_seg] int<lower=0>     seg_count;
  vector<lower=0>[n_seg]        seg_l;

  // ---- Spatial GP / HSGP block --------------------------------------------
  int<lower=1> S;                          // species (= 1 in this model)
  int<lower=1> M;                          // HSGP basis dimension = m1*m2*m3

  // Coordinates of segment midpoints, normalised to [-1, 1]
  // (col 1 = X km, col 2 = Y km, col 3 = Z_bathy m)
  matrix[n_seg, 3] coords;

  // Original-unit half-ranges used for normalisation
  vector<lower=0>[3] coord_scale;

  // HSGP boundary extension factors and basis counts
  vector<lower=0>[3]    L_hsgp;
  array[3] int<lower=1> m_hsgp;

  // Prior hyperparameters for the spatial field
  real prior_mu_sp_mu;
  real<lower=0> prior_mu_sp_sig;
  real<lower=0> prior_gp_sigma_shape;
  real<lower=0> prior_gp_sigma_rate;
  real prior_gp_lx_mu;       real<lower=0> prior_gp_lx_sig;
  real prior_gp_ly_mu;       real<lower=0> prior_gp_ly_sig;
  real prior_gp_lz_mu;       real<lower=0> prior_gp_lz_sig;

}

// =============================================================================
transformed data {

  // HSGP design matrix at segment midpoints — built once
  matrix[n_seg, M] PHI = hsgp_phi(coords, L_hsgp, m_hsgp);

  // Centred group sizes (per detection)
  vector[n] s_c = to_vector(s) - s_centre;

}

// =============================================================================
parameters {

  // Spatial field (S = 1)
  vector[S]              mu_sp;        // log group-density intercept
  vector<lower=0>[S]     gp_sigma;     // GP marginal SD per species
  matrix<lower=0>[S, 3]  gp_l;         // length-scales (lx km, ly km, lz m)
  matrix[S, M]           z_beta;       // non-centred basis coefficients

  // Detection function
  real log_sigma;                      // log(sigma) at s = s_centre
  real beta_size;                      // log_sigma slope on (s - s_centre)

  // Group-size sub-models (both always sampled; only one enters the
  // likelihood, the other sits at its prior — same convention as v4.1)
  real<lower=0> mu_s;
  real<lower=0> phi_s;
  real          mu_log_s;
  real<lower=0> sigma_log_s;

}

// =============================================================================
transformed parameters {

  // Per-species (S = 1) latent field at segments
  vector[n_seg]      f_seg;
  vector[n_seg]      log_lambda_groups;   // log group density at each segment
  vector<lower=0>[n_seg] lambda_groups;

  {
    vector[M] wt = hsgp_weights(
      gp_sigma[1],
      to_vector(gp_l[1, ]),
      L_hsgp,
      coord_scale,
      m_hsgp
    );
    f_seg = PHI * (wt .* to_vector(z_beta[1, ]));
    for (j in 1:n_seg) {
      // Clamp log_lambda for numerical safety (matches v3.2's [-10, 15])
      log_lambda_groups[j] = fmax(fmin(mu_sp[1] + f_seg[j], 15.0), -10.0);
      lambda_groups[j]     = exp(log_lambda_groups[j]);
    }
  }

  // Per-detection HN scale + ESW
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

  // Representative sigma / ESW at population mean group size (reporting
  // only; not used in encounter-rate likelihood — that uses esw_pop)
  real<lower=0> sigma_rep;
  if (use_size_covar == 1)
    sigma_rep = exp(log_sigma + beta_size * (mu_s - s_centre));
  else
    sigma_rep = exp(log_sigma);
  real<lower=0> esw_rep = sigma_rep * sqrt(pi() / 2)
                          * erf(w / (sqrt(2) * sigma_rep));

  // Population-corrected ESW. Identical to v4.1 — sums
  // sigma(k)-induced ESW(k) over the group-size pmf P(k) for k=1..S_max.
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
      real sigma_k = exp(log_sigma + beta_size * (k - s_centre));
      real esw_k   = sigma_k * sqrt(pi() / 2)
                     * erf(w / (sqrt(2) * sigma_k));
      esw_sum += pk * esw_k;
      p_sum   += pk;
    }
    esw_pop = esw_sum / p_sum;
  } else {
    esw_pop = esw_rep;
  }

}

// =============================================================================
model {

  // ----- Spatial-field priors -------------------------------------------
  mu_sp       ~ normal(prior_mu_sp_mu, prior_mu_sp_sig);
  gp_sigma    ~ gamma(prior_gp_sigma_shape, prior_gp_sigma_rate);
  for (sp in 1:S) {
    gp_l[sp, 1] ~ normal(prior_gp_lx_mu, prior_gp_lx_sig);
    gp_l[sp, 2] ~ normal(prior_gp_ly_mu, prior_gp_ly_sig);
    gp_l[sp, 3] ~ normal(prior_gp_lz_mu, prior_gp_lz_sig);
  }
  to_vector(z_beta) ~ std_normal();

  // ----- Detection-function + group-size priors -------------------------
  log_sigma    ~ normal(log_sigma_prior_mean,    log_sigma_prior_sd);
  beta_size    ~ normal(beta_size_prior_mean,    beta_size_prior_sd);
  mu_s         ~ gamma(mu_s_prior_shape,         mu_s_prior_rate);
  phi_s        ~ gamma(phi_s_prior_shape,        phi_s_prior_rate);
  mu_log_s     ~ normal(mu_log_prior_mean,       mu_log_prior_sd);
  sigma_log_s  ~ gamma(sigma_log_prior_shape,    sigma_log_prior_rate);

  // ----- Half-normal detection-distance likelihood (per detection) ------
  for (i in 1:n) {
    target += -square(x[i]) / (2.0 * square(sigma_i[i])) - log(esw_i[i]);
  }

  // ----- Group-size likelihood with size-bias correction -----------------
  if (group_size_dist == 0) {
    target += neg_binomial_2_lpmf(s | mu_s, phi_s)
              - n * log1m(neg_binomial_2_cdf(0 | mu_s, phi_s));
  } else {
    target += lognormal_lpdf(s | mu_log_s, sigma_log_s);
  }
  target += sum(log(esw_i)) - n * log(esw_pop);

  // ----- SPATIAL Poisson encounter-rate likelihood -----------------------
  // Per-segment expected detected groups:
  //   E[seg_count[j]] = strip_area_j * lambda_groups[j] * (esw_pop / w)
  //                   = (2 * w * seg_l[j]) * lambda_groups[j] * (esw_pop / w)
  //                   = 2 * seg_l[j] * esw_pop * lambda_groups[j]
  target += poisson_lpmf(seg_count | 2.0 * seg_l .* lambda_groups * esw_pop);

}

// =============================================================================
generated quantities {

  // Mean detection probability in the strip (population-corrected) and at mu_s
  real<lower=0, upper=1> p       = esw_pop / w;
  real<lower=0, upper=1> p_at_mu = esw_rep / w;

  // Spatial PPCs ---------------------------------------------------------
  // Per-segment fitted lambda + animal density (mu_s scales groups -> animals)
  vector[n_seg] lambda_animals_seg;
  for (j in 1:n_seg) {
    lambda_animals_seg[j] = lambda_groups[j] * mu_s;
  }

  // Posterior-predictive segment counts (PPC for the encounter-rate model)
  array[n_seg] int seg_count_rep;
  for (j in 1:n_seg) {
    real mu_j = 2.0 * seg_l[j] * esw_pop * lambda_groups[j];
    // Cap mu to avoid Stan poisson_rng overflow at extreme draws
    real mu_safe = fmin(mu_j, 1e7);
    seg_count_rep[j] = poisson_rng(mu_safe);
  }

  // Per-segment encounter-rate log-lik (for loo / PSIS-LOO)
  vector[n_seg] log_lik_seg;
  for (j in 1:n_seg) {
    real mu_j = 2.0 * seg_l[j] * esw_pop * lambda_groups[j];
    log_lik_seg[j] = poisson_lpmf(seg_count[j] | mu_j);
  }

  // Per-detection HN log-lik (for loo / PSIS-LOO on the distance side)
  vector[n] log_lik_dist;
  for (i in 1:n) {
    log_lik_dist[i] = -square(x[i]) / (2.0 * square(sigma_i[i]))
                      - log(esw_i[i]);
  }

  // Animal-density summary (mean of lambda_animals across segments)
  real D_animals_mean = mean(lambda_animals_seg);

}
