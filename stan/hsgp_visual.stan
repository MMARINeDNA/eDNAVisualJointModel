// =============================================================================
// hsgp_visual.stan
//
// Spatial line-transect (visual survey) model for ONE species:
//
//   log lambda(x_j) = mu_sp + f(X_j, Y_j) + B_bathy_j . beta_bathy      animals / km^2
//   f ~ GP(0, K_SE(lx, ly)),  2-D HSGP over normalised (X, Y)
//   seg_count_j ~ Poisson(lambda_j / E[s] * 2 * L_j * esw_pop)
//   x_i ~ half-normal(sigma_i) truncated at w,  log sigma_i = log_sigma + beta_size (s_i - s_centre)
//   s_i ~ group-size distribution, size-bias corrected for detection
//
// The latent field has the same structure and units as the eDNA model
// (HSGP4eDNA stan/hsgp_2d_bathysp.stan), so the two can share it in the
// joint model. The detection / group-size machinery is the validated
// distance_hn_dens_v4.1.stan, moved into stan/include/visual_functions.stan.
//
// Switches (data):
//   use_gp  = 0  drops the spatial field (f = 0; GP parameters have size 0)
//   K_bathy = 0  drops the bottom-depth spline
// With both off this is the non-spatial v4.1 model, re-parameterised on animal
// density (lambda = exp(mu_sp)) instead of group density.
//
// Compile with
//   cmdstan_model("stan/hsgp_visual.stan",
//                 include_paths = c("stan/include", "external/HSGP4eDNA/stan/include"))
// =============================================================================

functions {
#include hsgp_functions.stan
#include visual_functions.stan
}

data {
  // ---- Detections (one row per detected group) ------------------------------
  int<lower=1> n;
  vector<lower=0>[n] x;                    // perpendicular distances (km)
  array[n] int<lower=1> s;                 // group sizes
  real<lower=0> w;                         // truncation distance (km)

  // ---- Detection function ---------------------------------------------------
  real log_sigma_prior_mean;
  real<lower=0> log_sigma_prior_sd;
  int<lower=0, upper=1> use_size_covar;
  real s_centre;
  real beta_size_prior_mean;
  real<lower=0> beta_size_prior_sd;

  // ---- Group-size distribution ---------------------------------------------
  int<lower=1> S_max;                      // pmf support 1..S_max
  int<lower=0, upper=1> group_size_dist;   // 0 = ZT neg-binomial, 1 = log-normal
  real<lower=0> mu_s_prior_shape;
  real<lower=0> mu_s_prior_rate;
  real<lower=0> phi_s_prior_shape;
  real<lower=0> phi_s_prior_rate;
  real mu_log_prior_mean;
  real<lower=0> mu_log_prior_sd;
  real<lower=0> sigma_log_prior_shape;
  real<lower=0> sigma_log_prior_rate;

  // ---- Segments (encounter rate) --------------------------------------------
  int<lower=1> n_seg;
  array[n_seg] int<lower=0> seg_count;     // detected groups per segment
  vector<lower=0>[n_seg] seg_l;            // segment length (km)

  // ---- Latent field ---------------------------------------------------------
  int<lower=0, upper=1> use_gp;
  int<lower=1> D1;                         // GP dimensions (2: X, Y)
  int<lower=1> M;                          // number of basis functions
  array[M, D1] int INDICES;                // basis index per axis
  matrix[n_seg, D1] coords;                // segment midpoints, normalised to [-1, 1]
  vector<lower=0>[D1] coord_scale;         // half-range per axis (km)
  array[D1] real<lower=0> L_hsgp;          // boundary factor per axis
  int<lower=0> K_bathy;                    // bottom-depth spline basis size
  matrix[n_seg, K_bathy] B_bathy;

  int<lower=0> N_pred;                     // prediction locations
  matrix[N_pred, D1] pred_coords;
  matrix[N_pred, K_bathy] B_bathy_pred;

  // ---- Field priors ---------------------------------------------------------
  real prior_mu_sp_mu;
  real<lower=0> prior_mu_sp_sig;
  real<lower=0> prior_gp_sigma_shape;
  real<lower=0> prior_gp_sigma_rate;
  // gp_l_raw[d] (normalised units: ell / coord_scale) ~ inv_gamma(shape, scale),
  // per axis; set by R/gp_priors.R gp_prior_from_design()
  vector<lower=0>[D1] prior_gp_l_shape;
  vector<lower=0>[D1] prior_gp_l_scale;
  real<lower=0> prior_beta_bathy_sig;
}

transformed data {
  matrix[n_seg, M]  PHI      = hsgp_phi(L_hsgp, INDICES, coords);
  matrix[N_pred, M] PHI_pred = hsgp_phi(L_hsgp, INDICES, pred_coords);
  vector[S_max] k_support    = linspaced_vector(S_max, 1, S_max);
}

parameters {
  // Latent field
  real mu_sp;                              // log animal density intercept
  vector<lower=0>[use_gp] gp_sigma;
  matrix<lower=0>[use_gp, D1] gp_l_raw;    // length-scales, normalised units
  matrix[use_gp, M] z_beta;                // non-centred basis coefficients
  vector[K_bathy] beta_bathy;

  // Detection function
  real log_sigma;                          // log(sigma) at s = s_centre
  real beta_size;

  // Group size (both families sampled; only group_size_dist's enters the data
  // likelihood, the other sits at its prior - same convention as v4.1)
  real<lower=0> mu_s;
  real<lower=0> phi_s;
  real mu_log_s;
  real<lower=0> sigma_log_s;
}

transformed parameters {
  vector[n_seg] f_seg = rep_vector(0, n_seg);
  if (use_gp == 1) {
    f_seg = PHI * (hsgp_sqrt_spd(gp_sigma[1], gp_l_raw[1], L_hsgp, INDICES, D1)
                   .* z_beta[1]');
  }
  vector[n_seg] log_lambda = mu_sp + f_seg + B_bathy * beta_bathy;

  vector[S_max] pk = group_size_pmf(S_max, group_size_dist, mu_s, phi_s,
                                    mu_log_s, sigma_log_s);
  real<lower=0> esw_pop = esw_population(pk, use_size_covar, log_sigma,
                                         beta_size, s_centre, w);
  real<lower=1> mean_group_size = dot_product(k_support, pk);
}

model {
  // Field priors
  mu_sp ~ normal(prior_mu_sp_mu, prior_mu_sp_sig);
  gp_sigma ~ gamma(prior_gp_sigma_shape, prior_gp_sigma_rate);
  for (g in 1:use_gp) {
    gp_l_raw[g]' ~ inv_gamma(prior_gp_l_shape, prior_gp_l_scale);
  }
  to_vector(z_beta) ~ std_normal();
  beta_bathy ~ normal(0, prior_beta_bathy_sig);

  // Detection + group-size priors (as v4.1)
  log_sigma   ~ normal(log_sigma_prior_mean, log_sigma_prior_sd);
  beta_size   ~ normal(beta_size_prior_mean, beta_size_prior_sd);
  mu_s        ~ gamma(mu_s_prior_shape, mu_s_prior_rate);
  phi_s       ~ gamma(phi_s_prior_shape, phi_s_prior_rate);
  mu_log_s    ~ normal(mu_log_prior_mean, mu_log_prior_sd);
  sigma_log_s ~ gamma(sigma_log_prior_shape, sigma_log_prior_rate);

  // Likelihood
  target += sum(lt_detection_loglik(x, s, w, use_size_covar, log_sigma,
                                    beta_size, s_centre, group_size_dist,
                                    mu_s, phi_s, mu_log_s, sigma_log_s, esw_pop));
  target += sum(lt_segment_loglik(seg_count, seg_l, log_lambda,
                                  log(mean_group_size), esw_pop));
}

generated quantities {
  row_vector[D1] gp_l = rep_row_vector(0, D1);   // length-scales (km); 0 if use_gp = 0
  if (use_gp == 1) gp_l = gp_l_raw[1] .* coord_scale';

  real p_det    = esw_pop / w;                     // mean detection prob. in strip
  real sigma_c  = hn_sigma(log_sigma);             // HN scale at s = s_centre (km)
  real D_mean   = mean(exp(log_lambda));           // mean animal density over segments

  vector[n]     log_lik_det = lt_detection_loglik(x, s, w, use_size_covar, log_sigma,
                                beta_size, s_centre, group_size_dist,
                                mu_s, phi_s, mu_log_s, sigma_log_s, esw_pop);
  vector[n_seg] log_lik_seg = lt_segment_loglik(seg_count, seg_l, log_lambda,
                                log(mean_group_size), esw_pop);

  array[n_seg] int pp_seg_count;
  for (j in 1:n_seg) {
    pp_seg_count[j] = poisson_log_rng(log_lambda[j] - log(mean_group_size)
                                      + log(2 * seg_l[j] * esw_pop));
  }

  vector[N_pred] log_lambda_pred = mu_sp + B_bathy_pred * beta_bathy;
  if (use_gp == 1) {
    log_lambda_pred += PHI_pred * (hsgp_sqrt_spd(gp_sigma[1], gp_l_raw[1], L_hsgp,
                                                 INDICES, D1) .* z_beta[1]');
  }
}
