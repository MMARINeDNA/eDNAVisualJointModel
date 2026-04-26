// =============================================================================
// whale_edna_hsgp.stan  (v3.2)
//
// 3-D anisotropic HSGP model — single species (Pacific hake), surface
// samples only, qPCR observations only. Stripped-down debug case for
// the v3 pipeline.
//
// Latent field:
//   log(lambda_i) = mu + f(X_i, Y_i, Z_bathy_i)   // HSGP GP residual
//
// Observation model — qPCR hurdle (hake only):
//   detect[r]  ~ Bernoulli(1 - exp(-kappa * lambda_edna[i] * vol_aliquot/100))
//   Ct[r]      ~ Normal(alpha_ct - beta_ct * log(lambda_edna[i] * vol_aliquot/100),
//                       sigma_ct)   if detect[r] == 1
//
// Differences from v3:
//   * S = 1 (single species). Vector / matrix shapes still indexed by S so
//     the structure mirrors v3, but in practice S is hard-coded to 1 by
//     the upstream R pipeline.
//   * Metabarcoding likelihood and all MB-only parameters
//     (beta0_phi / gamma0_phi / gamma1_phi) are removed.
//   * `kappa` is fixed in the data block (was sampled in v3). The
//     calibration is treated as known, mirroring alpha_ct / beta_ct.
//   * Numerical guards from v3 carried over verbatim:
//       - log_lam clamped to [-10, 15] before exp / linear use
//       - sigma_ct floored at 1e-6
//       - p_det squashed to [1e-9, 1-1e-9]
// =============================================================================

functions {

  // -------------------------------------------------------------------------
  // 1-D Laplacian eigenfunction on [-L, L]
  // -------------------------------------------------------------------------
  real phi1d(real L, int m, real x) {
    return (1.0 / sqrt(L)) * sin(m * pi() * (x + L) / (2.0 * L));
  }

  // -------------------------------------------------------------------------
  // Spectral density of 1-D squared-exponential kernel
  // -------------------------------------------------------------------------
  real spd_se_1d(real alpha, real rho, real w) {
    return square(alpha) * sqrt(2.0 * pi()) * rho
             * exp(-0.5 * square(rho * w));
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

  // Dimensions
  int<lower=1> N;          // total samples (one per station, surface only)
  int<lower=1> S;          // species (= 1 in v3.2)
  int<lower=1> M;          // HSGP basis dimension = m1 * m2 * m3

  // Conversion factor (log scale): log(copies per animal per litre per km^2)
  real log_conv_factor;

  // Coordinates normalised to [-1, 1] (col 1=X km, 2=Y km, 3=Z_bathy m)
  matrix[N, 3] coords;

  // Original-unit half-ranges used for normalisation
  vector<lower=0>[3] coord_scale;

  // HSGP boundary extension factors and basis counts
  vector<lower=0>[3]    L_hsgp;
  array[3] int<lower=1> m_hsgp;

  // Fixed water-column eDNA log-offsets (pre-computed in R): N x S.
  // For v3.2 (surface only) this is all zeros, kept for shape compatibility.
  matrix[N, S] log_zsample_effect;

  // Aliquot volume (microlitres); dilution fraction = vol_aliquot / 100
  real<lower=0> vol_aliquot;

  // ------ qPCR data (hake only, s=1), long form ------
  int<lower=1>                          N_qpcr_long;
  array[N_qpcr_long] int<lower=1, upper=N> qpcr_sample_idx;
  array[N_qpcr_long] int<lower=0, upper=1> qpcr_detect;
  vector[N_qpcr_long]                      qpcr_ct;      // 0 if not detected

  // Fixed qPCR standard-curve coefficients (pre-estimated, e.g. from the
  // hake survey standard curve). Treated as known.
  real alpha_ct;
  real beta_ct;

  // Fixed qPCR detection rate. Promoted from a parameter to data in
  // v3.2 (the calibration is treated as known, same as alpha_ct/beta_ct).
  real<lower=0, upper=1> kappa;

  // Prior hyperparameters (passed as data for easy tuning).
  // gp_sigma now uses a Gamma(shape, rate) prior. The previous
  // half_normal(0, prior_gp_sigma_sig) had its mode at 0, which - even
  // with informative qPCR Ct data - pinned the gp_sigma posterior near
  // zero (chains visibly trapped at the boundary, posterior mean three
  // orders of magnitude below truth). Gamma puts zero density at 0, so
  // the prior cannot trap the chain there; with shape=4, rate=2 the
  // mode sits at 1.5 and the mean at 2.0, comfortably bracketing the
  // simulated truth of 1.2.
  real prior_mu_sp_mu;
  real prior_mu_sp_sig;
  real<lower=0> prior_gp_sigma_shape;
  real<lower=0> prior_gp_sigma_rate;
  real prior_gp_lx_mu;       real<lower=0> prior_gp_lx_sig;
  real prior_gp_ly_mu;       real<lower=0> prior_gp_ly_sig;
  real prior_gp_lz_mu;       real<lower=0> prior_gp_lz_sig;
  real prior_sigma_ct_mu;    real<lower=0> prior_sigma_ct_sig;

}

// =============================================================================
transformed data {

  // HSGP design matrix - built once, reused every iteration
  matrix[N, M] PHI = hsgp_phi(coords, L_hsgp, m_hsgp);

  // Log aliquot dilution: log(vol_aliquot / 100)
  real log_vol_frac = log(vol_aliquot / 100.0);

}

// =============================================================================
parameters {

  vector[S]             mu_sp;      // log-density intercepts
  vector<lower=0>[S]    gp_sigma;   // GP marginal SD per species
  matrix<lower=0>[S, 3] gp_l;       // GP length-scales: lx(km), ly(km), lz(m)
  matrix[S, M]          z_beta;     // non-centred basis coefficients

  // qPCR Ct noise. alpha_ct, beta_ct, kappa are fixed in the data block.
  real<lower=0> sigma_ct;

}

// =============================================================================
transformed parameters {

  matrix[N, S] log_lambda;       // log true animal density
  matrix[N, S] log_lambda_edna;  // log eDNA at sample depth

  for (s in 1:S) {
    vector[M] wt  = hsgp_weights(
      gp_sigma[s],
      to_vector(gp_l[s, ]),
      L_hsgp,
      coord_scale,
      m_hsgp
    );
    vector[N] f_s = PHI * (wt .* to_vector(z_beta[s, ]));

    for (i in 1:N) {
      log_lambda[i, s]      = mu_sp[s] + f_s[i];
      log_lambda_edna[i, s] = log_lambda[i, s]
                              + log_zsample_effect[i, s]
                              + log_conv_factor;
    }
  }

}

// =============================================================================
model {

  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------
  mu_sp    ~ normal(prior_mu_sp_mu, prior_mu_sp_sig);
  gp_sigma ~ gamma(prior_gp_sigma_shape, prior_gp_sigma_rate);

  for (s in 1:S) {
    gp_l[s, 1] ~ normal(prior_gp_lx_mu, prior_gp_lx_sig);
    gp_l[s, 2] ~ normal(prior_gp_ly_mu, prior_gp_ly_sig);
    gp_l[s, 3] ~ normal(prior_gp_lz_mu, prior_gp_lz_sig);
  }

  to_vector(z_beta) ~ std_normal();

  sigma_ct ~ normal(prior_sigma_ct_mu, prior_sigma_ct_sig);

  // ------------------------------------------------------------------
  // Likelihood - qPCR hurdle (hake only, s=1)
  // ------------------------------------------------------------------
  {
    for (r in 1:N_qpcr_long) {
      int  i            = qpcr_sample_idx[r];
      // log copies in aliquot = log eDNA at depth + log dilution fraction
      real log_lam      = log_lambda_edna[i, 1] + log_vol_frac;
      real lam          = exp(fmin(log_lam, 15.0));   // guard overflow

      real p_det        = 1.0 - exp(-kappa * lam);
      p_det             = fmax(fmin(p_det, 1.0 - 1e-9), 1e-9);

      target += bernoulli_lpmf(qpcr_detect[r] | p_det);

      if (qpcr_detect[r] == 1) {
        real log_lam_safe = fmax(fmin(log_lam, 15.0), -10.0);
        real mu_ct        = alpha_ct - beta_ct * log_lam_safe;
        target += normal_lpdf(qpcr_ct[r] | mu_ct, fmax(sigma_ct, 1e-6));
      }
    }
  }

}

// =============================================================================
generated quantities {

  vector[N_qpcr_long] log_lik_qpcr;

  array[N_qpcr_long] int  pp_qpcr_detect;
  vector[N_qpcr_long]     pp_qpcr_ct;

  matrix[N, S] lambda_hat;
  matrix[N, S] lambda_edna_hat;

  for (i in 1:N) {
    for (s in 1:S) {
      lambda_hat[i, s]      = exp(log_lambda[i, s]);
      lambda_edna_hat[i, s] = exp(log_lambda_edna[i, s]);
    }
  }

  for (r in 1:N_qpcr_long) {
    int  i            = qpcr_sample_idx[r];
    real log_lam      = log_lambda_edna[i, 1] + log_vol_frac;
    real lam          = exp(fmin(log_lam, 15.0));
    real p_det        = fmax(fmin(1.0 - exp(-kappa * lam), 1.0 - 1e-9), 1e-9);

    log_lik_qpcr[r] = bernoulli_lpmf(qpcr_detect[r] | p_det);
    if (qpcr_detect[r] == 1) {
      real log_lam_safe = fmax(fmin(log_lam, 15.0), -10.0);
      real mu_ct        = alpha_ct - beta_ct * log_lam_safe;
      log_lik_qpcr[r]  += normal_lpdf(qpcr_ct[r] | mu_ct, fmax(sigma_ct, 1e-6));
    }

    pp_qpcr_detect[r] = bernoulli_rng(p_det);
    if (pp_qpcr_detect[r] == 1) {
      real log_lam_safe = fmax(fmin(log_lam, 15.0), -10.0);
      real mu_ct        = alpha_ct - beta_ct * log_lam_safe;
      pp_qpcr_ct[r]     = normal_rng(mu_ct, fmax(sigma_ct, 1e-6));
    } else {
      pp_qpcr_ct[r] = 0.0;
    }
  }

}
