// =============================================================================
// hsgp_visual.stan
//
// Spatial line-transect (visual survey) model for ONE species:
//
//   log lambda(x_j) = mu_sp + f(X_j, Y_j) + B_bathy_j . beta_bathy      animals / km^2
//   f ~ GP(0, K_SE(lx, ly)),  2-D HSGP over normalised (X, Y)
//   seg_count_j ~ Poisson(lambda_j / E[s] * 2 * L_j * esw_pop(b_j))
//                 (or negative binomial with size phi_nb, use_nb = 1)
//   x_i ~ half-normal or hazard-rate (det_key) detection, scale sigma_i, truncated at w,
//         log sigma_i = log_sigma + beta_size (s_i - s_centre) + beta_bft (b_i - b_centre)
//   s_i ~ group-size distribution, size-bias corrected for detection
// b = Beaufort sea state of the segment (binned to levels; bft_level holds
// b - b_centre per level), so ESW differs between segments.
//
// The latent field has the same structure and units as the eDNA model
// (HSGP4eDNA stan/hsgp_2d_bathysp.stan), so the two can share it in the
// joint model. The detection / group-size machinery is the validated
// distance_hn_dens_v4.1.stan, moved into stan/include/visual_functions.stan.
//
// Switches (data):
//   use_gp  = 0  drops the spatial field (f = 0; GP parameters have size 0)
//   K_bathy = 0  drops the bottom-depth spline
//   use_bft_covar = 0  drops the sea-state term (n_bft = 1, bft_level = 0):
//                      exactly the single-ESW model
//   use_nb  = 1  negative-binomial segment counts (clustered sightings), with
//                1 / sqrt(phi_nb) ~ half-normal(0, prior_nb_inv_sqrt_sd);
//                use_nb = 0 is exactly the Poisson model (phi_nb has size 0)
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
  int<lower=0, upper=1> det_key;           // 0 = half-normal, 1 = hazard-rate
  real log_b_prior_mean;                   // hazard-rate shape b, log scale
  real<lower=0> log_b_prior_sd;
  // Sea-state (Beaufort) covariate on log sigma. Segments are binned into
  // n_bft levels; the size-averaged ESW is computed once per level.
  int<lower=0, upper=1> use_bft_covar;
  int<lower=1> n_bft;
  vector[n_bft] bft_level;                 // Beaufort - centre, per level
  real beta_bft_prior_mean;
  real<lower=0> beta_bft_prior_sd;

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
  array[n_seg] int<lower=1, upper=n_bft> seg_bft;  // sea-state level per segment
  array[n] int<lower=1, upper=n_bft> det_bft;      // ... per detection (its segment's)
  int<lower=0, upper=1> use_nb;            // 1 = negative-binomial segment counts
  real<lower=0> prior_nb_inv_sqrt_sd;      // half-normal sd on 1 / sqrt(phi_nb)

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
  int Q_hr = 400;                          // midpoint cells for the hazard-rate ESW
  vector[Q_hr] log_xq = log((linspaced_vector(Q_hr, 1, Q_hr) - 0.5) * w / Q_hr);
  // Group-size bins for the ESW sum: sizes 1..30 individually, then 10% log bins
  array[S_max] int ksz_bin;
  for (k in 1:S_max) {
    ksz_bin[k] = k <= 30 ? k : 31 + to_int(floor(log(k) / log(1.1))) - to_int(floor(log(31) / log(1.1)));
  }
  int n_kbin = ksz_bin[S_max];
}

parameters {
  // Latent field
  real mu_sp;                              // log animal density intercept
  vector<lower=0>[use_gp] gp_sigma;
  matrix<lower=0>[use_gp, D1] gp_l_raw;    // length-scales, normalised units
  matrix[use_gp, M] z_beta;                // non-centred basis coefficients
  vector[K_bathy] beta_bathy;

  // Detection function
  real log_sigma;                          // log(sigma) at s = s_centre, Beaufort = centre
  real beta_size;
  vector[use_bft_covar] beta_bft;          // change in log sigma per Beaufort unit
  vector[det_key] log_b;                   // hazard-rate shape (size 0 if half-normal)

  // Segment-count overdispersion (size 0 unless use_nb)
  vector<lower=0>[use_nb] inv_sqrt_phi_nb;

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
  real b_hr = 1.0;
  if (det_key == 1) b_hr = exp(log_b[1]);
  vector[n_bft] log_sigma_bft = rep_vector(log_sigma, n_bft);
  if (use_bft_covar == 1) log_sigma_bft += beta_bft[1] * bft_level;
  vector<lower=0>[n_bft] esw_bft;          // population-average ESW per sea-state level
  for (b in 1:n_bft) {
    esw_bft[b] = esw_population(pk, use_size_covar, det_key, log_sigma_bft[b],
                                beta_size, s_centre, b_hr, w, log_xq, ksz_bin, n_kbin);
  }
  real<lower=1> mean_group_size = dot_product(k_support, pk);
  real phi_nb = 1e300;                      // NB size; unused unless use_nb
  if (use_nb == 1) phi_nb = inv_square(inv_sqrt_phi_nb[1]);
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
  beta_bft    ~ normal(beta_bft_prior_mean, beta_bft_prior_sd);
  log_b       ~ normal(log_b_prior_mean, log_b_prior_sd);
  inv_sqrt_phi_nb ~ normal(0, prior_nb_inv_sqrt_sd);
  mu_s        ~ gamma(mu_s_prior_shape, mu_s_prior_rate);
  phi_s       ~ gamma(phi_s_prior_shape, phi_s_prior_rate);
  mu_log_s    ~ normal(mu_log_prior_mean, mu_log_prior_sd);
  sigma_log_s ~ gamma(sigma_log_prior_shape, sigma_log_prior_rate);

  // Likelihood
  target += sum(lt_detection_loglik(x, s, w, use_size_covar, det_key, log_sigma_bft[det_bft],
                                    beta_size, s_centre, b_hr, group_size_dist,
                                    mu_s, phi_s, mu_log_s, sigma_log_s, esw_bft[det_bft]));
  if (use_nb == 1) {
    target += sum(lt_segment_loglik_nb(seg_count, seg_l, log_lambda,
                                       log(mean_group_size), esw_bft[seg_bft], phi_nb));
  } else {
    target += sum(lt_segment_loglik(seg_count, seg_l, log_lambda,
                                    log(mean_group_size), esw_bft[seg_bft]));
  }
}

generated quantities {
  row_vector[D1] gp_l = rep_row_vector(0, D1);   // length-scales (km); 0 if use_gp = 0
  if (use_gp == 1) gp_l = gp_l_raw[1] .* coord_scale';

  real p_det    = mean(esw_bft[seg_bft]) / w;      // mean detection prob. in strip (over segments)
  real sigma_c  = hn_sigma(log_sigma);             // key scale at s = s_centre, Beaufort centre (km)
  real b_shape  = b_hr;                            // hazard-rate shape (1, unused, if half-normal)
  real D_mean   = mean(exp(log_lambda));           // mean animal density over segments

  vector[n]     log_lik_det = lt_detection_loglik(x, s, w, use_size_covar, det_key,
                                log_sigma_bft[det_bft], beta_size, s_centre, b_hr, group_size_dist,
                                mu_s, phi_s, mu_log_s, sigma_log_s, esw_bft[det_bft]);
  vector[n_seg] log_lik_seg = use_nb == 1
    ? lt_segment_loglik_nb(seg_count, seg_l, log_lambda, log(mean_group_size), esw_bft[seg_bft], phi_nb)
    : lt_segment_loglik(seg_count, seg_l, log_lambda, log(mean_group_size), esw_bft[seg_bft]);

  array[n_seg] int pp_seg_count;
  for (j in 1:n_seg) {
    real eta_j = log_lambda[j] - log(mean_group_size) + log(2 * seg_l[j] * esw_bft[seg_bft[j]]);
    pp_seg_count[j] = use_nb == 1 ? neg_binomial_2_log_rng(eta_j, phi_nb) : poisson_log_rng(eta_j);
  }

  vector[N_pred] log_lambda_pred = mu_sp + B_bathy_pred * beta_bathy;
  if (use_gp == 1) {
    log_lambda_pred += PHI_pred * (hsgp_sqrt_spd(gp_sigma[1], gp_l_raw[1], L_hsgp,
                                                 INDICES, D1) .* z_beta[1]');
  }
}
