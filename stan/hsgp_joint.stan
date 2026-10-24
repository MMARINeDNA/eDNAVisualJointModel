// =============================================================================
// hsgp_joint.stan
//
// Joint eDNA + visual line-transect model. ONE latent animal-density field per
// species, shared by both observation models:
//
//   log lambda_s(x) = mu_s + f_s(X, Y) + B(Z_bathy) . beta_s       animals / km^2
//   f_s ~ GP(0, K_SE(lx_s, ly_s)), 2-D HSGP over (X, Y) normalised by the DOMAIN
//
// eDNA (all S species + junk; HSGP4eDNA v1.12 gp2d_bathysp's observation model,
// without the mechanistic phi):
//   log lambda_edna_is = log lambda_s(x_i) + log_zsample_effect_is
//                        + B_zdepth_i . beta_zdepth_s   (sampling-depth smooth)
//                        + log_conv_s + log_vol_bottle_i + log g_is
//   -> qPCR hurdle (species 1) + zero-inflated beta-binomial metabarcoding,
//      g_is between-bottle noise (one shared theta), g_i,junk junk noise
//      (own theta_junk). Each qPCR well and metabarcoding replicate adds its
//      log dilution to the aliquot fraction. The qPCR standard curve is fixed
//      (qpcr_curve_model = 0) or, per plate, estimated jointly from the
//      standards and the field wells (qpcr_curve_model = 1).
//   With no dilutions, every bottle at log_vol_filtered, the fixed curve,
//   ct_sd_cond = 0 and K_zdepth = 0 this is the eDNA model of HSGP4eDNA v1.5
//   exactly (tests/joint_exactness.R).
//
// Visual (the L line-transect species, lt_species_idx maps them into 1..S;
// the hsgp_visual.stan observation model):
//   seg_count_lj ~ Poisson(lambda_l(x_j) / E[s_l] * 2 L_j * esw_l(b_j)), or
//                  negative binomial with size phi_nb_l for species with
//                  use_nb[l] = 1 (clustered sightings); use_nb[l] = 0 is the
//                  Poisson model exactly
//   detections: half-normal or hazard-rate distances (det_key; optional
//   group-size covariate) + size-bias-corrected group sizes, with
//   log sigma = log_sigma_l + beta_size_l (s - s_centre) + beta_bft_l (b - b_centre)
//               + beta_vis_l (v - v_centre)
//   b = Beaufort sea state and v = visibility of the segment; segments are
//   grouped into n_bft detection-condition cells shared by all species
//   (bft_level, vis_level = cell value - centre). use_bft_covar[l] = 0 /
//   use_vis_covar[l] = 0 drops the term for species l, giving its model
//   without it exactly.
//
// Switches (data):
//   use_edna = 0    drop the eDNA likelihood (visual-only fit)
//   use_visual = 0  drop the visual likelihood (eDNA-only fit)
//   n_conv_free > 0 estimate log_conv for the listed species (with the
//                   visual data pinning absolute density) instead of fixing it
// The same model, priors and basis serve all three fits, so eDNA-only vs
// visual-only vs joint comparisons differ only in the data used.
//
// Compile with
//   cmdstan_model("stan/hsgp_joint.stan",
//                 include_paths = c("stan/include", "external/HSGP4eDNA/stan/include"))
// =============================================================================

functions {
#include hsgp_functions.stan
#include edna_functions.stan
#include visual_functions.stan
}

data {
  // ---- Switches ---------------------------------------------------------------
  int<lower=0, upper=1> use_edna;
  int<lower=0, upper=1> use_visual;

  // ---- Latent field (shared) --------------------------------------------------
  int<lower=1> S;                          // species
  int<lower=1> D1;                         // GP dimensions (2: X, Y)
  int<lower=1> M;                          // basis functions
  array[M, D1] int INDICES;
  array[D1] real<lower=0> L_hsgp;
  vector<lower=0>[D1] coord_scale;         // km per normalised unit
  int<lower=0> K_bathy;                    // bottom-depth spline basis size

  real prior_gp_sigma_shape;
  real<lower=0> prior_gp_sigma_rate;
  vector<lower=0>[D1] prior_gp_l_shape;    // inv_gamma on gp_l_raw (design-based)
  vector<lower=0>[D1] prior_gp_l_scale;
  vector[S + 1] prior_mu_sp_mu;            // + junk
  vector<lower=0>[S + 1] prior_mu_sp_sig;
  real<lower=0> prior_beta_bathy_sig;

  // ---- eDNA ---------------------------------------------------------------------
  int<lower=1> N;                          // samples (station x depth)
  matrix[N, D1] coords;                    // normalised by the domain
  matrix[N, K_bathy] B_bathy;
  matrix[N, S] log_zsample_effect;
  vector[S + 1] log_conv_factor;           // fixed values (+ junk)
  real log_vol_filtered;
  real<lower=0> vol_aliquot;

  int<lower=0> N_qpcr_long;
  array[N_qpcr_long] int<lower=1, upper=N> qpcr_sample_idx;
  array[N_qpcr_long] int<lower=0, upper=1> qpcr_detect;
  vector[N_qpcr_long] qpcr_ct;
  int<lower=0> N_mb_long;
  array[N_mb_long] int<lower=1, upper=N> mb_sample_idx;
  array[N_mb_long, S + 1] int<lower=0> mb_reads;
  array[N_mb_long] int<lower=0> mb_total;

  real alpha_ct;
  real beta_ct;
  real<lower=0, upper=1> kappa;
  real gamma0_ct;
  real gamma1_ct;
  real<lower=0> sigma0_ct;

  real prior_beta0_phi_mu;   real<lower=0> prior_beta0_phi_sig;
  real prior_gamma0_phi_mu;  real<lower=0> prior_gamma0_phi_sig;
  real prior_gamma1_phi_mu;  real<lower=0> prior_gamma1_phi_sig;

  // Per-well qPCR dilution, per-replicate metabarcoding dilution and filtered
  // volume per bottle (HSGP4eDNA v1.8): log offsets on the aliquot copies
  vector[N_qpcr_long] qpcr_log_dilution;  // log fraction of template (0 = undiluted)
  vector[N_mb_long] mb_log_dilution;      // -log fold (0 = undiluted)
  vector[N] log_vol_bottle;               // log litres filtered per bottle

  // qPCR standard curves estimated in the model (HSGP4eDNA v1.9).
  // qpcr_curve_model = 0: the fixed curve above. 1: intercept and slope per
  // plate (slope >= beta_ct_min), Ct-noise terms and detection rate shared,
  // from the standards and the field wells.
  int<lower=0, upper=1> qpcr_curve_model;
  int<lower=1> N_plate;
  array[N_qpcr_long] int<lower=1, upper=N_plate> qpcr_plate_idx;
  int<lower=0> N_std;
  array[N_std] int<lower=1, upper=N_plate> std_plate_idx;
  array[N_std] int<lower=0, upper=1> std_detect;
  vector[N_std] std_ct;                   // 0 if not detected
  vector[N_std] std_log_copies;           // log(template uL x copies per uL)
  vector[2] prior_alpha_ct;
  vector[2] prior_beta_ct;
  real<lower=0> beta_ct_min;              // 1 / log(2)
  vector[2] prior_gamma0_ct;
  vector[2] prior_gamma1_ct;
  real<lower=0> prior_sigma0_ct_sig;
  vector[2] prior_logit_kappa;
  // Ct-noise SD at very low copies (HSGP4eDNA v1.12): 0 = at the expected
  // copies; 1 = at the expected copies given detection, bounded by one copy
  int<lower=0, upper=1> ct_sd_cond;

  // Sampling depth as a per-species smooth on the eDNA copies (HSGP4eDNA
  // v1.10; K_zdepth = 0: none)
  int<lower=0> K_zdepth;
  matrix[N, K_zdepth] B_zdepth;
  real<lower=0> prior_beta_zdepth_sig;

  // Between-bottle and junk noise (HSGP4eDNA v1.5 semantics)
  int<lower=0, upper=1> use_bottle_noise;
  int<lower=1, upper=2> bottle_noise_dist;   // 1 gamma (exact NB), 2 lognormal
  int<lower=0, upper=1> bottle_theta_shared;
  real<lower=0> prior_cv_bottle_sig;
  int<lower=0, upper=1> use_junk_noise;
  int<lower=1, upper=2> junk_noise_dist;     // 1 gamma (centred), 2 lognormal (non-centred)
  real<lower=0> prior_cv_junk_sig;

  // Free conversion factors (science step): species whose log_conv is estimated
  int<lower=0, upper=S> n_conv_free;
  array[n_conv_free] int<lower=1, upper=S> conv_free_idx;
  vector[n_conv_free] prior_log_conv_mu;
  vector<lower=0>[n_conv_free] prior_log_conv_sig;

  // ---- Visual line transects ------------------------------------------------------
  int<lower=0> L;                          // line-transect species
  array[L] int<lower=1, upper=S> lt_species_idx;
  int<lower=0> n_seg;
  matrix[n_seg, D1] seg_coords;            // normalised by the domain
  matrix[n_seg, K_bathy] B_bathy_seg;
  vector<lower=0>[n_seg] seg_l;            // km
  array[L, n_seg] int<lower=0> seg_count;  // detected groups per segment
  array[L] real<lower=0> w;                // truncation distance (km)

  array[L] int<lower=0> n_det;             // detections per species ...
  int<lower=0> n_det_total;                // ... stacked species by species
  vector<lower=0>[n_det_total] det_x;      // perpendicular distance (km)
  array[n_det_total] int<lower=1> det_s;   // group size

  vector[L] log_sigma_prior_mean;
  vector<lower=0>[L] log_sigma_prior_sd;
  array[L] int<lower=0, upper=1> use_size_covar;
  vector[L] s_centre;
  vector[L] beta_size_prior_mean;
  vector<lower=0>[L] beta_size_prior_sd;
  array[L] int<lower=0, upper=1> det_key;  // 0 = half-normal, 1 = hazard-rate
  vector[L] log_b_prior_mean;              // hazard-rate shape b, log scale
  vector<lower=0>[L] log_b_prior_sd;
  // Sea-state (Beaufort) covariate on log sigma, per species. Segments are
  // binned into n_bft levels shared by all species; the size-averaged ESW is
  // computed once per species and level.
  array[L] int<lower=0, upper=1> use_bft_covar;
  int<lower=1> n_bft;
  vector[n_bft] bft_level;                 // Beaufort - centre, per level
  array[n_seg] int<lower=1, upper=n_bft> seg_bft;        // level per segment
  array[n_det_total] int<lower=1, upper=n_bft> det_bft;  // ... per detection (its segment's)
  vector[L] beta_bft_prior_mean;
  vector<lower=0>[L] beta_bft_prior_sd;
  // Visibility covariate on log sigma, per species, over the same cells
  array[L] int<lower=0, upper=1> use_vis_covar;
  vector[n_bft] vis_level;                 // visibility - centre, per cell
  vector[L] beta_vis_prior_mean;
  vector<lower=0>[L] beta_vis_prior_sd;
  // Negative-binomial segment counts per species; 1 / sqrt(phi_nb) ~
  // half-normal(0, prior_nb_inv_sqrt_sd)
  array[L] int<lower=0, upper=1> use_nb;
  real<lower=0> prior_nb_inv_sqrt_sd;
  array[L] int<lower=1> S_max;
  array[L] int<lower=0, upper=1> group_size_dist;
  vector<lower=0>[L] mu_s_prior_shape;
  vector<lower=0>[L] mu_s_prior_rate;
  vector<lower=0>[L] phi_s_prior_shape;
  vector<lower=0>[L] phi_s_prior_rate;
  vector[L] mu_log_prior_mean;
  vector<lower=0>[L] mu_log_prior_sd;
  vector<lower=0>[L] sigma_log_prior_shape;
  vector<lower=0>[L] sigma_log_prior_rate;
}

transformed data {
  matrix[N, M]     PHI     = hsgp_phi(L_hsgp, INDICES, coords);
  matrix[n_seg, M] PHI_seg = hsgp_phi(L_hsgp, INDICES, seg_coords);
  int n_cv = use_bottle_noise * (bottle_theta_shared == 1 ? 1 : S);
  real log_vol_frac = log(vol_aliquot / 100.0);
  // Per-replicate log offsets on the aliquot copies: the aliquot's fraction of
  // the extract plus its dilution
  vector[N_qpcr_long] qpcr_log_offset = log_vol_frac + qpcr_log_dilution;
  vector[N_mb_long] mb_log_offset = log_vol_frac + mb_log_dilution;
  array[L] int det_start;
  {
    int pos = 1;
    for (l in 1:L) { det_start[l] = pos; pos += n_det[l]; }
  }
  // Detection-function quadrature and group-size bins, as hsgp_visual.stan:
  // midpoint cells for the hazard-rate ESW, and sizes 1..30 individually then
  // 10% log bins for the size-averaged ESW sum.
  int Q_hr = 400;
  array[L] vector[Q_hr] log_xq;
  int S_max_all = L > 0 ? max(S_max) : 1;
  array[L, S_max_all] int ksz_bin = rep_array(1, L, S_max_all);
  array[L] int n_kbin;
  for (l in 1:L) {
    log_xq[l] = log((linspaced_vector(Q_hr, 1, Q_hr) - 0.5) * w[l] / Q_hr);
    for (k in 1:S_max[l])
      ksz_bin[l, k] = k <= 30 ? k : 31 + to_int(floor(log(k) / log(1.1))) - to_int(floor(log(31) / log(1.1)));
    n_kbin[l] = ksz_bin[l, S_max[l]];
  }
  // Hazard-rate species get a shape parameter
  int n_hr = sum(det_key);
  array[L] int hr_pos = rep_array(0, L);
  {
    int j = 0;
    for (l in 1:L) if (det_key[l] == 1) { j += 1; hr_pos[l] = j; }
  }
  // Species with the sea-state covariate get a Beaufort slope
  int n_bft_sp = sum(use_bft_covar);
  array[L] int bft_pos = rep_array(0, L);
  {
    int j = 0;
    for (l in 1:L) if (use_bft_covar[l] == 1) { j += 1; bft_pos[l] = j; }
  }
  // Species with the visibility covariate get a visibility slope
  int n_vis_sp = sum(use_vis_covar);
  array[L] int vis_pos = rep_array(0, L);
  {
    int j = 0;
    for (l in 1:L) if (use_vis_covar[l] == 1) { j += 1; vis_pos[l] = j; }
  }
  // Species with negative-binomial counts get an overdispersion parameter
  int n_nb_sp = sum(use_nb);
  array[L] int nb_pos = rep_array(0, L);
  {
    int j = 0;
    for (l in 1:L) if (use_nb[l] == 1) { j += 1; nb_pos[l] = j; }
  }
}

parameters {
  // Latent field
  vector[S + 1] mu_sp;                     // + junk intercept
  vector<lower=0>[S] gp_sigma;
  matrix<lower=0>[S, D1] gp_l_raw;
  array[S] vector[M] z_beta;
  matrix[S, K_bathy] beta_bathy;

  // eDNA observation model
  real          beta0_phi;
  real<lower=0> gamma0_phi;
  real<lower=0> gamma1_phi;
  matrix[N * use_bottle_noise, S] g_raw;
  vector<lower=0>[n_cv] cv_bottle;
  vector[N * use_junk_noise] g_junk_raw;
  vector<lower=0>[use_junk_noise] cv_junk;
  vector[n_conv_free] log_conv_free;
  // qPCR standard curves (qpcr_curve_model = 1)
  vector[N_plate * qpcr_curve_model] alpha_ct_plate;
  vector<lower=beta_ct_min>[N_plate * qpcr_curve_model] beta_ct_plate;
  vector[qpcr_curve_model] gamma0_ct_est;
  vector[qpcr_curve_model] gamma1_ct_est;
  vector<lower=0>[qpcr_curve_model] sigma0_ct_est;
  vector[qpcr_curve_model] logit_kappa_est;
  matrix[S, K_zdepth] beta_zdepth;         // sampling-depth spline coefficients

  // Visual observation model, per line-transect species
  vector[L] log_sigma;
  vector[L] beta_size;
  vector[n_hr] log_b;                      // hazard-rate shape (hazard-rate species only)
  vector[n_bft_sp] beta_bft;               // change in log sigma per Beaufort unit (bft species only)
  vector[n_vis_sp] beta_vis;               // change in log sigma per visibility unit (vis species only)
  vector<lower=0>[n_nb_sp] inv_sqrt_phi_nb; // segment-count overdispersion (NB species only)
  vector<lower=0>[L] mu_s;
  vector<lower=0>[L] phi_s;
  vector[L] mu_log_s;
  vector<lower=0>[L] sigma_log_s;
}

transformed parameters {
  matrix[N, S]     log_lambda;               // field at eDNA samples
  matrix[n_seg, S] log_lambda_seg;           // field at transect segments
  for (s in 1:S) {
    vector[M] w_s = hsgp_sqrt_spd(gp_sigma[s], gp_l_raw[s, ], L_hsgp, INDICES, D1) .* z_beta[s];
    log_lambda[, s]     = mu_sp[s] + PHI * w_s     + B_bathy     * beta_bathy[s]';
    log_lambda_seg[, s] = mu_sp[s] + PHI_seg * w_s + B_bathy_seg * beta_bathy[s]';
  }

  vector[S + 1] log_conv = log_conv_factor;
  log_conv[conv_free_idx] = log_conv_free;

  matrix[N, S + 1] log_lambda_edna;
  for (s in 1:S) {
    log_lambda_edna[, s] = log_lambda[, s] + log_zsample_effect[, s]
                           + log_conv[s] + log_vol_bottle;
    if (K_zdepth > 0) log_lambda_edna[, s] += B_zdepth * beta_zdepth[s]';
    if (use_bottle_noise == 1)
      log_lambda_edna[, s] += bottle_dist_log_g(g_raw[, s],
                                cv_bottle[bottle_theta_shared == 1 ? 1 : s], bottle_noise_dist);
  }
  log_lambda_edna[, S + 1] = mu_sp[S + 1] + log_conv[S + 1] + log_vol_bottle;
  if (use_junk_noise == 1)
    log_lambda_edna[, S + 1] += bottle_dist_log_g(g_junk_raw, cv_junk[1], junk_noise_dist);

  // Visual: population-average ESW per species and sea-state level, and mean
  // group size per species. Without the covariate every level shares one ESW.
  vector[L] b_hr = rep_vector(1.0, L);     // hazard-rate shape (1, unused, if half-normal)
  for (l in 1:L) if (det_key[l] == 1) b_hr[l] = exp(log_b[hr_pos[l]]);
  array[L] vector[n_bft] log_sigma_bft;    // log sigma at s = s_centre, per level
  array[L] vector[n_bft] esw_bft;
  vector[L] mean_group_size;
  for (l in 1:L) {
    vector[S_max[l]] pk = group_size_pmf(S_max[l], group_size_dist[l], mu_s[l], phi_s[l],
                                         mu_log_s[l], sigma_log_s[l]);
    log_sigma_bft[l] = rep_vector(log_sigma[l], n_bft);
    if (use_bft_covar[l] == 1 || use_vis_covar[l] == 1) {
      if (use_bft_covar[l] == 1) log_sigma_bft[l] += beta_bft[bft_pos[l]] * bft_level;
      if (use_vis_covar[l] == 1) log_sigma_bft[l] += beta_vis[vis_pos[l]] * vis_level;
      for (b in 1:n_bft)
        esw_bft[l][b] = esw_population(pk, use_size_covar[l], det_key[l], log_sigma_bft[l][b],
                                       beta_size[l], s_centre[l], b_hr[l], w[l], log_xq[l],
                                       ksz_bin[l, 1:S_max[l]], n_kbin[l]);
    } else {
      esw_bft[l] = rep_vector(esw_population(pk, use_size_covar[l], det_key[l], log_sigma[l],
                                             beta_size[l], s_centre[l], b_hr[l], w[l], log_xq[l],
                                             ksz_bin[l, 1:S_max[l]], n_kbin[l]), n_bft);
    }
    mean_group_size[l] = dot_product(linspaced_vector(S_max[l], 1, S_max[l]), pk);
  }
}

model {
  // ---- Field priors ----
  mu_sp ~ normal(prior_mu_sp_mu, prior_mu_sp_sig);
  gp_sigma ~ gamma(prior_gp_sigma_shape, prior_gp_sigma_rate);
  for (s in 1:S) {
    gp_l_raw[s, ]' ~ inv_gamma(prior_gp_l_shape, prior_gp_l_scale);
    z_beta[s] ~ std_normal();
  }
  to_vector(beta_bathy) ~ normal(0, prior_beta_bathy_sig);

  // ---- eDNA priors ----
  beta0_phi  ~ normal(prior_beta0_phi_mu,  prior_beta0_phi_sig);
  gamma0_phi ~ normal(prior_gamma0_phi_mu, prior_gamma0_phi_sig);
  gamma1_phi ~ normal(prior_gamma1_phi_mu, prior_gamma1_phi_sig);
  cv_bottle ~ normal(0, prior_cv_bottle_sig);
  for (s in 1:(S * use_bottle_noise)) {
    if (bottle_noise_dist == 1)
      g_raw[, s] ~ log_gamma_mixing(1 / square(cv_bottle[bottle_theta_shared == 1 ? 1 : s]));
    else
      g_raw[, s] ~ std_normal();
  }
  cv_junk ~ normal(0, prior_cv_junk_sig);
  if (use_junk_noise == 1) {
    if (junk_noise_dist == 1) g_junk_raw ~ log_gamma_mixing(1 / square(cv_junk[1]));
    else                        g_junk_raw ~ std_normal();
  }
  log_conv_free ~ normal(prior_log_conv_mu, prior_log_conv_sig);
  alpha_ct_plate  ~ normal(prior_alpha_ct[1], prior_alpha_ct[2]);
  beta_ct_plate   ~ normal(prior_beta_ct[1], prior_beta_ct[2]);
  gamma0_ct_est   ~ normal(prior_gamma0_ct[1], prior_gamma0_ct[2]);
  gamma1_ct_est   ~ normal(prior_gamma1_ct[1], prior_gamma1_ct[2]);
  sigma0_ct_est   ~ normal(0, prior_sigma0_ct_sig);
  logit_kappa_est ~ normal(prior_logit_kappa[1], prior_logit_kappa[2]);
  to_vector(beta_zdepth) ~ normal(0, prior_beta_zdepth_sig);

  // ---- Visual priors ----
  log_sigma   ~ normal(log_sigma_prior_mean, log_sigma_prior_sd);
  beta_size   ~ normal(beta_size_prior_mean, beta_size_prior_sd);
  for (l in 1:L) if (det_key[l] == 1) log_b[hr_pos[l]] ~ normal(log_b_prior_mean[l], log_b_prior_sd[l]);
  for (l in 1:L) if (use_bft_covar[l] == 1)
    beta_bft[bft_pos[l]] ~ normal(beta_bft_prior_mean[l], beta_bft_prior_sd[l]);
  for (l in 1:L) if (use_vis_covar[l] == 1)
    beta_vis[vis_pos[l]] ~ normal(beta_vis_prior_mean[l], beta_vis_prior_sd[l]);
  inv_sqrt_phi_nb ~ normal(0, prior_nb_inv_sqrt_sd);
  mu_s        ~ gamma(mu_s_prior_shape, mu_s_prior_rate);
  phi_s       ~ gamma(phi_s_prior_shape, phi_s_prior_rate);
  mu_log_s    ~ normal(mu_log_prior_mean, mu_log_prior_sd);
  sigma_log_s ~ gamma(sigma_log_prior_shape, sigma_log_prior_rate);

  // ---- eDNA likelihood ----
  // (the standards inform the curves only when the eDNA data are used)
  if (use_edna == 1) {
    if (qpcr_curve_model == 0) {
      target += sum(edna_qpcr_loglik_offset(qpcr_sample_idx, qpcr_detect, qpcr_ct,
                                            log_lambda_edna[, 1], qpcr_log_offset,
                                            alpha_ct, beta_ct, kappa,
                                            gamma0_ct, gamma1_ct, sigma0_ct, ct_sd_cond));
    } else {
      // field wells and standards, through the same per-plate curves
      real kappa_q = inv_logit(logit_kappa_est[1]);
      target += sum(edna_qpcr_plate_loglik(qpcr_detect, qpcr_ct,
                      log_lambda_edna[qpcr_sample_idx, 1] + qpcr_log_offset, qpcr_plate_idx,
                      alpha_ct_plate, beta_ct_plate, kappa_q,
                      gamma0_ct_est[1], gamma1_ct_est[1], sigma0_ct_est[1], ct_sd_cond));
      target += sum(edna_qpcr_plate_loglik(std_detect, std_ct, std_log_copies, std_plate_idx,
                      alpha_ct_plate, beta_ct_plate, kappa_q,
                      gamma0_ct_est[1], gamma1_ct_est[1], sigma0_ct_est[1], ct_sd_cond));
    }
    target += sum(edna_mb_loglik_offset(mb_sample_idx, mb_reads, mb_total,
                                        log_lambda_edna, mb_log_offset,
                                        beta0_phi, gamma0_phi, gamma1_phi));
  }

  // ---- Visual likelihood ----
  if (use_visual == 1) {
    for (l in 1:L) {
      int n_l = n_det[l];
      if (n_l > 0) {
        array[n_l] int db = segment(det_bft, det_start[l], n_l);
        target += sum(lt_detection_loglik(
          segment(det_x, det_start[l], n_l), segment(det_s, det_start[l], n_l), w[l],
          use_size_covar[l], det_key[l], log_sigma_bft[l][db], beta_size[l],
          s_centre[l], b_hr[l], group_size_dist[l], mu_s[l], phi_s[l], mu_log_s[l], sigma_log_s[l],
          esw_bft[l][db]));
      }
      if (use_nb[l] == 1) {
        target += sum(lt_segment_loglik_nb(seg_count[l], seg_l,
                                           log_lambda_seg[, lt_species_idx[l]],
                                           log(mean_group_size[l]),
                                           esw_bft[l][seg_bft],
                                           inv_square(inv_sqrt_phi_nb[nb_pos[l]])));
      } else {
        target += sum(lt_segment_loglik(seg_count[l], seg_l,
                                        log_lambda_seg[, lt_species_idx[l]],
                                        log(mean_group_size[l]),
                                        esw_bft[l][seg_bft]));
      }
    }
  }
}

generated quantities {
  matrix[S, D1] gp_l;                      // length-scales (km)
  for (d in 1:D1) gp_l[, d] = gp_l_raw[, d] * coord_scale[d];
  vector[n_cv] theta_bottle = bottle_noise_dist == 1
      ? 1 ./ square(cv_bottle) : 1 ./ (exp(square(cv_bottle)) - 1);
  vector[use_junk_noise] theta_junk = junk_noise_dist == 1
      ? 1 ./ square(cv_junk) : 1 ./ (exp(square(cv_junk)) - 1);
  vector[L] sigma_c = exp(log_sigma);      // key scale at s = s_centre, Beaufort centre (km)
  vector[L] b_shape = b_hr;                // hazard-rate shape (1 if half-normal)
  vector[L] p_det;                         // mean detection prob. in the strip (over segments)
  for (l in 1:L) p_det[l] = mean(esw_bft[l][seg_bft]) / w[l];
  vector[n_nb_sp] phi_nb = inv_square(inv_sqrt_phi_nb);   // NB size (NB species only)
  // qPCR curves in the eDNA model's terms (qpcr_curve_model = 1)
  vector[qpcr_curve_model] kappa_qpcr = inv_logit(logit_kappa_est);
  vector[N_plate * qpcr_curve_model] qpcr_efficiency = exp(1 ./ beta_ct_plate) - 1;  // 1 = doubling
}
