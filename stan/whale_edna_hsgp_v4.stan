// =============================================================================
// whale_edna_hsgp.stan  (v4)
//
// 3-D anisotropic HSGP model for three-species whale/fish eDNA.
// SF to US/Canada border, UTM Zone 10N.
//
// Latent field:
//   log(lambda_si) = mu_s
//                  + f_s(X_i, Y_i, Z_bathy_i)   // HSGP GP residual
//
// Observation models (unchanged from v3) — qPCR hurdle on Pacific hake
// and zero-inflated Beta-Binomial metabarcoding on all three species.
// alpha_ct / beta_ct are fixed data; kappa / sigma_ct / phi-hyperparameters
// retain Normal priors with mu/sigma passed as data.
//
// What changed relative to v3:
//   * No structural change to the Stan model — it already consumed data
//     in long form (one row per aliquot) indexed by qpcr_sample_idx /
//     mb_sample_idx. That representation naturally accommodates the v4
//     simulation's variable per-sample replication (most samples have 3
//     reps but ~20% have 1 or 2), so N_qpcr_long and N_mb_long are no
//     longer assumed to equal N * 3.
//
// Fixes carried from v3:
//   1. log1m_exp argument clamped to keep it strictly negative
//   2. alpha_bb / beta_bb clamped to >= 1e-6
//   3. More conservative init values (handled in runner)
//   4. lam_sum clamped before log(lam_sum) to avoid log(0)
// =============================================================================

functions {

  // -------------------------------------------------------------------------
  // Zero-inflated Beta-Binomial log-PMF
  // log_p_zero = log P(absent from bottle) = -lambda  (Poisson)
  // log_p_pos  = log(1 - exp(log_p_zero))
  // count == 0 : log_sum_exp(log_p_zero, log_p_pos + BB(0|n,a,b))
  // count  > 0 : log_p_pos + BB(count|n,a,b)
  // -------------------------------------------------------------------------
  real zi_beta_binomial_lpmf(int count, int reads,
                              real log_p_zero, real log_p_pos,
                              real alpha_bb, real beta_bb) {
    real log_bb = beta_binomial_lpmf(count | reads, alpha_bb, beta_bb);
    if (count == 0) {
      return log_sum_exp(log_p_zero, log_p_pos + log_bb);
    } else {
      return log_p_pos + log_bb;
    }
  }

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
  int<lower=1> N;          // total samples (station × depth combinations)
  int<lower=1> S;          // species (3)
  int<lower=1> M;          // HSGP basis dimension = m1 * m2 * m3

  // Species-specific conversion factor (log scale): log(copies per animal
  // per litre per km^2). Whales shed more eDNA per animal than hake, so
  // this is per-species.
  vector[S] log_conv_factor;

  // Coordinates normalised to [-1, 1] (col 1=X km, 2=Y km, 3=Z_bathy m)
  matrix[N, 3] coords;

  // Original-unit half-ranges used for normalisation
  vector<lower=0>[3] coord_scale;

  // HSGP boundary extension factors and basis counts
  vector<lower=0>[3]    L_hsgp;
  array[3] int<lower=1> m_hsgp;

  // Fixed water-column eDNA log-offsets (pre-computed in R): N × S
  matrix[N, S] log_zsample_effect;

  // Aliquot volume (microlitres); dilution fraction = vol_aliquot / 100
  real<lower=0> vol_aliquot;

  // ------ qPCR data (hake only, s=1), long form ------
  int<lower=1>                          N_qpcr_long;
  array[N_qpcr_long] int<lower=1, upper=N> qpcr_sample_idx;
  array[N_qpcr_long] int<lower=0, upper=1> qpcr_detect;
  vector[N_qpcr_long]                      qpcr_ct;      // 0 if not detected

  // ------ Metabarcoding data (all species), long form ------
  int<lower=1>                          N_mb_long;
  array[N_mb_long] int<lower=1, upper=N> mb_sample_idx;
  array[N_mb_long, S] int<lower=0>       mb_reads;
  array[N_mb_long]    int<lower=0>       mb_total;

  // 1 = zero-inflated BB, 0 = plain BB
  int<lower=0, upper=1> use_zi;

  // Fixed qPCR standard-curve coefficients (pre-estimated, e.g. from the
  // hake survey standard curve). These were parameters in earlier versions
  // but are now passed as data.
  real alpha_ct;
  real beta_ct;

  // Prior hyperparameters (passed as data for easy tuning)
  real prior_mu_sp_mu;
  real prior_mu_sp_sig;
  real prior_gp_sigma_sig;
  real prior_gp_lx_mu;       real prior_gp_lx_sig;
  real prior_gp_ly_mu;       real prior_gp_ly_sig;
  real prior_gp_lz_mu;       real prior_gp_lz_sig;
  real prior_kappa_mu;       real<lower=0> prior_kappa_sig;
  real prior_sigma_ct_mu;    real<lower=0> prior_sigma_ct_sig;
  real prior_beta0_phi_mu;   real<lower=0> prior_beta0_phi_sig;
  real prior_gamma0_phi_mu;  real<lower=0> prior_gamma0_phi_sig;
  real prior_gamma1_phi_mu;  real<lower=0> prior_gamma1_phi_sig;

}

// =============================================================================
transformed data {

  // HSGP design matrix — built once, reused every iteration
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

  // qPCR (hake only). alpha_ct and beta_ct are fixed (see data block).
  real<lower=0, upper=1>   kappa;
  real<lower=0>            sigma_ct;

  // Metabarcoding overdispersion
  vector[S]          beta0_phi;
  vector<lower=0>[S] gamma0_phi;
  vector<lower=0>[S] gamma1_phi;

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
                              + log_conv_factor[s];
    }
  }

}

// =============================================================================
model {

  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------
  mu_sp    ~ normal(prior_mu_sp_mu, prior_mu_sp_sig);
  gp_sigma ~ normal(0, prior_gp_sigma_sig);

  for (s in 1:S) {
    gp_l[s, 1] ~ normal(prior_gp_lx_mu, prior_gp_lx_sig);
    gp_l[s, 2] ~ normal(prior_gp_ly_mu, prior_gp_ly_sig);
    gp_l[s, 3] ~ normal(prior_gp_lz_mu, prior_gp_lz_sig);
  }

  to_vector(z_beta) ~ std_normal();

  kappa    ~ normal(prior_kappa_mu,    prior_kappa_sig);
  sigma_ct ~ normal(prior_sigma_ct_mu, prior_sigma_ct_sig);

  beta0_phi  ~ normal(prior_beta0_phi_mu,  prior_beta0_phi_sig);
  gamma0_phi ~ normal(prior_gamma0_phi_mu, prior_gamma0_phi_sig);
  gamma1_phi ~ normal(prior_gamma1_phi_mu, prior_gamma1_phi_sig);

  // ------------------------------------------------------------------
  // Likelihood 1 — qPCR hurdle (hake only, s=1)
  // ------------------------------------------------------------------
  {
    for (r in 1:N_qpcr_long) {
      int  i            = qpcr_sample_idx[r];
      // log copies in aliquot = log eDNA at depth + log dilution fraction
      real log_lam      = log_lambda_edna[i, 1] + log_vol_frac;
      real lam          = exp(fmin(log_lam, 15.0));   // FIX 3: guard overflow

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

  // ------------------------------------------------------------------
  // Likelihood 2 — Metabarcoding ZI-Beta-Binomial (all species)
  // ------------------------------------------------------------------
  {
    for (r in 1:N_mb_long) {
      int i = mb_sample_idx[r];

      // eDNA concentrations at this sample
      vector[S] lam_edna_i;
      for (s in 1:S) lam_edna_i[s] = exp(fmin(log_lambda_edna[i, s], 15.0));

      // FIX 4: clamp lam_sum before taking log to avoid log(0)
      real lam_sum      = fmax(sum(lam_edna_i), 1e-12);
      real log_lam_sum  = log(lam_sum);

      // Compositional proportions
      vector[S] pi_i;
      if (lam_sum < 1e-12) {
        pi_i = rep_vector(1.0 / S, S);
      } else {
        pi_i = lam_edna_i / lam_sum;
      }

      for (s in 1:S) {
        // FIX 1: clamp log_lam_s before exp to keep log_p_zero well-defined
        real log_lam_s  = fmax(fmin(log_lambda_edna[i, s], 15.0), -15.0);
        real lam_s      = exp(log_lam_s);

        real log_phi_s  = beta0_phi[s] + log_lam_sum
                          + fmax(gamma0_phi[s] - gamma1_phi[s] * log_lam_s, 0.0);
        real phi_s      = exp(fmin(log_phi_s, 10.0));   // cap phi to avoid overflow

        real pi_s       = fmax(fmin(pi_i[s], 1.0 - 1e-6), 1e-6);

        // FIX 2: clamp alpha_bb and beta_bb away from zero
        real alpha_bb   = fmax(pi_s * phi_s,         1e-6);
        real beta_bb    = fmax((1.0 - pi_s) * phi_s, 1e-6);

        if (use_zi) {
          // FIX 1 (cont): keep log1m_exp argument strictly negative
          real log_p_zero = -lam_s;
          real log_p_pos  = log1m_exp(fmin(log_p_zero, -1e-9));

          target += zi_beta_binomial_lpmf(
            mb_reads[r, s] | mb_total[r],
            log_p_zero, log_p_pos,
            alpha_bb, beta_bb
          );
        } else {
          target += beta_binomial_lpmf(
            mb_reads[r, s] | mb_total[r],
            alpha_bb, beta_bb
          );
        }
      }
    }
  }

}

// =============================================================================
generated quantities {

  vector[N_qpcr_long] log_lik_qpcr;
  vector[N_mb_long]   log_lik_mb;

  array[N_qpcr_long] int  pp_qpcr_detect;
  vector[N_qpcr_long]     pp_qpcr_ct;
  array[N_mb_long, S] int pp_mb_reads;

  matrix[N, S] lambda_hat;
  matrix[N, S] lambda_edna_hat;

  for (i in 1:N) {
    for (s in 1:S) {
      lambda_hat[i, s]      = exp(log_lambda[i, s]);
      lambda_edna_hat[i, s] = exp(log_lambda_edna[i, s]);
    }
  }

  // ------------------------------------------------------------------
  // qPCR generated quantities
  // ------------------------------------------------------------------
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

  // ------------------------------------------------------------------
  // Metabarcoding generated quantities
  // ------------------------------------------------------------------
  for (r in 1:N_mb_long) {
    log_lik_mb[r] = 0.0;

    if (mb_total[r] == 0) {
      for (s in 1:S) pp_mb_reads[r, s] = 0;
      continue;
    }

    int i = mb_sample_idx[r];

    vector[S] lam_edna_i;
    for (s in 1:S) lam_edna_i[s] = exp(fmin(log_lambda_edna[i, s], 15.0));

    // FIX 4: clamp before log
    real lam_sum     = fmax(sum(lam_edna_i), 1e-12);
    real log_lam_sum = log(lam_sum);

    vector[S] pi_i;
    if (lam_sum < 1e-12) {
      pi_i = rep_vector(1.0 / S, S);
    } else {
      pi_i = lam_edna_i / lam_sum;
    }

    for (s in 1:S) {
      real log_lam_s  = fmax(fmin(log_lambda_edna[i, s], 15.0), -15.0);
      real lam_s      = exp(log_lam_s);

      real log_phi_s  = beta0_phi[s] + log_lam_sum
                        + fmax(gamma0_phi[s] - gamma1_phi[s] * log_lam_s, 0.0);
      real phi_s      = exp(fmin(log_phi_s, 10.0));

      real pi_s       = fmax(fmin(pi_i[s], 1.0 - 1e-6), 1e-6);

      // FIX 2: clamp BB parameters
      real alpha_bb   = fmax(pi_s * phi_s,         1e-6);
      real beta_bb    = fmax((1.0 - pi_s) * phi_s, 1e-6);

      if (use_zi) {
        // FIX 1: keep log1m_exp argument strictly negative
        real log_p_zero = -lam_s;
        real log_p_pos  = log1m_exp(fmin(log_p_zero, -1e-9));

        log_lik_mb[r] += zi_beta_binomial_lpmf(
          mb_reads[r, s] | mb_total[r],
          log_p_zero, log_p_pos,
          alpha_bb, beta_bb
        );
      } else {
        log_lik_mb[r] += beta_binomial_lpmf(
          mb_reads[r, s] | mb_total[r],
          alpha_bb, beta_bb
        );
      }

      real p_bb         = beta_rng(alpha_bb, beta_bb);
      pp_mb_reads[r, s] = binomial_rng(mb_total[r],
                            fmax(fmin(p_bb, 1.0 - 1e-9), 1e-9));
    }
  }

}
