// =============================================================================
// whale_edna_hsgp.stan  (v4.2)
//
// 3-D anisotropic HSGP model for three-species whale/fish eDNA.
// SF to US/Canada border, UTM Zone 10N.
//
// Latent field:
//   log(lambda_si)      = mu_s + f_s(X_i, Y_i, Z_bathy_i)   // HSGP GP
//   log(lambda_edna_si) = log_lambda_si
//                       + log_zsample_effect_si
//                       + log_conv_factor[s]
//                       + log_vol_filtered                  // <-- v4.1 fix
//
// Observation models - qPCR hurdle on Pacific hake (s = 1) and
// zero-inflated Beta-Binomial metabarcoding on all three species.
//
// Changes from v4.1 (applying v3.2 lessons):
//   * kappa, gamma_0, gamma_1, sigma_0 moved from `parameters` to `data`. The qPCR
//     standard-curve calibration is treated as known. sigma_ct is
//     a function of these parameters
//   * gp_sigma prior switched from half_normal(0, sig) to
//     gamma(shape, rate) - half_normal had its mode at 0 and was
//     visibly trapping gp_sigma posterior near zero in v3.2.
//   * log_vol_filtered added to log_lambda_edna - v3 / v4 model was
//     missing this factor, biasing mu_sp upward by log(vol_filtered).
//
// Fixes carried from v3 / v4:
//   1. log1m_exp argument clamped to keep it strictly negative
//   2. alpha_bb / beta_bb clamped to >= 1e-6
//   3. lam_sum clamped before log(lam_sum) to avoid log(0)
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

  // log(vol_filtered) - litres of seawater filtered. v4.1 fix: this
  // factor was missing from v3 / v4's log_lambda_edna, forcing mu_sp
  // upward by log(vol_filtered) to compensate. Now in the eDNA mean.
  real log_vol_filtered;

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
  array[N_mb_long, S+1] int<lower=0>       mb_reads; // Added +1 for junk
  array[N_mb_long]    int<lower=0>       mb_total;

  // 1 = zero-inflated BB, 0 = plain BB
  int<lower=0, upper=1> use_zi;

  // Fixed qPCR standard-curve coefficients (pre-estimated calibration).
  // alpha_ct, beta_ct were already data in v4. 
  // v4.2: kappa, gamma0_ct, gamma1_ct, sigma0_ct join them as data - the entire standard-curve fit is treated as
  // known. 
  real alpha_ct;
  real beta_ct;
  real<lower=0, upper=1> kappa;
  real gamma0_ct;
  real gamma1_ct;
  real<lower=0> sigma0_ct;

  // Prior hyperparameters (passed as data for easy tuning).
  // gp_sigma uses Gamma(shape, rate) - mode (alpha-1)/beta, mean
  // alpha/beta. With shape=8, rate=4: mode 1.75, mean 2.0, sd 0.71;
  // P(sigma < 0.3) ~ 4e-5 so the low-sigma trap from v3.2 is
  // unreachable.
  real prior_mu_sp_mu;
  real prior_mu_sp_sig;
  real<lower=0> prior_gp_sigma_shape;
  real<lower=0> prior_gp_sigma_rate;
  real prior_gp_lx_mu;       real<lower=0> prior_gp_lx_sig;
  real prior_gp_ly_mu;       real<lower=0> prior_gp_ly_sig;
  real prior_gp_lz_mu;       real<lower=0> prior_gp_lz_sig;
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

  vector[S+1]             mu_sp;      // log-density intercepts +1 = junk
  vector<lower=0>[S]    gp_sigma;   // GP marginal SD per species
  matrix<lower=0>[S, 3] gp_l;       // GP length-scales: lx(km), ly(km), lz(m)
  matrix[S, M]          z_beta;     // non-centred basis coefficients

  real<lower=0>         sigma_junk; // standard deviation among junk
  vector[N]             log_RE_junk_raw; // realizations of the junk RE among sites.
 
  // alpha_ct, beta_ct, kappa, sigma_ct are all fixed in the data block.

  // Metabarcoding overdispersion
  real          beta0_phi;
  real<lower=0> gamma0_phi;
  real<lower=0> gamma1_phi;

}

// =============================================================================
transformed parameters {

  matrix[N, S] log_lambda;       // log true animal density
  matrix[N, S+1] log_lambda_edna;  // log eDNA at sample depth +1=junk
  vector[N] log_RE_junk;  //  

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
                              + log_conv_factor[s]
                              + log_vol_filtered;
    }
  }
  
  // add latent variable for the junk for log_lambda_edna
  // NON-SPATIAL
  log_RE_junk = log_RE_junk_raw * sigma_junk;
  for (i in 1:N) {
      int J = S+1;
      log_lambda_edna[i, J] = mu_sp[J]
                              + log_RE_junk[i] 
                              + log_vol_filtered;
  }
}

// =============================================================================
model {

  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------
  mu_sp    ~ normal(prior_mu_sp_mu, prior_mu_sp_sig); // scalar mean and SD
  gp_sigma ~ gamma(prior_gp_sigma_shape, prior_gp_sigma_rate);

  log_RE_junk ~ std_normal(); // 
  sigma_junk ~ normal(0,2); // 
  
  for (s in 1:S) {
    gp_l[s, 1] ~ normal(prior_gp_lx_mu, prior_gp_lx_sig);
    gp_l[s, 2] ~ normal(prior_gp_ly_mu, prior_gp_ly_sig);
    gp_l[s, 3] ~ normal(prior_gp_lz_mu, prior_gp_lz_sig);
  }

  to_vector(z_beta) ~ std_normal();

  // alpha_ct, beta_ct, kappa, gamma0_ct, gamma1_ct, sigma0_ct are fixed in the data block (no prior).

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
        real sigma_ct     = pow(pow(sigma0_ct,2) + exp(2*(gamma0_ct + gamma1_ct*log_lam_safe )),0.5) ;
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
      vector[S+1] lam_edna_i;
      for (s in 1:(S+1)) lam_edna_i[s] = exp(fmin(log_lambda_edna[i, s], 15.0));

      // FIX 4: clamp lam_sum before taking log to avoid log(0)
      real lam_sum      = fmax(sum(lam_edna_i), 1e-12);
      real log_lam_sum  = log(lam_sum);

      // Compositional proportions
      vector[S+1] pi_i;
      // if (lam_sum < 1e-12) {
      //   pi_i = rep_vector(1.0 / (S+1), S+1);
      // } else {
        pi_i = lam_edna_i / lam_sum;
      // }

      for (s in 1:(S+1)) {
        // FIX 1: clamp log_lam_s before exp to keep log_p_zero well-defined
        real log_lam_s  = fmax(fmin(log_lambda_edna[i, s], 15.0), -15.0);
        real lam_s      = exp(log_lam_s);

        real log_phi_s  = beta0_phi + log_lam_sum
                          + fmax(gamma0_phi - gamma1_phi * log_lam_s, 0.0);
        real phi_s      = exp(fmin(log_phi_s, 10.0));   // cap phi to avoid overflow

        real pi_s       = fmax(fmin(pi_i[s], 1.0 - 1e-6), 1e-6);

        // FIX 2: clamp alpha_bb and beta_bb away from zero
        real alpha_bb   = fmax(pi_s * phi_s,         1e-6);
        real beta_bb    = fmax((1.0 - pi_s) * phi_s, 1e-6);

        if (use_zi) {
          // FIX 1 (cont): keep log1m_exp argument strictly negative
          real log_p_zero = -lam_s;
          real log_p_pos  = log1m_exp(fmin(log_p_zero, -1e-9));

          // target += zi_beta_binomial_lpmf(
          //   mb_reads[r, s] | mb_total[r],
          //   log_p_zero, log_p_pos,
          //   alpha_bb, beta_bb
          // );
        } else {
          // target += beta_binomial_lpmf(
          //   mb_reads[r, s] | mb_total[r],
          //   alpha_bb, beta_bb
          // );
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
  array[N_mb_long, S+1] int pp_mb_reads;

  matrix[N, S] lambda_hat;
  matrix[N, S+1] lambda_edna_hat;

  for (i in 1:N) {
    for (s in 1:S) {
      lambda_hat[i, s]      = exp(log_lambda[i, s]);
    }
    for (s in 1:(S+1)) {
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
      real sigma_ct     = pow(pow(sigma0_ct,2) + exp(2*(gamma0_ct + gamma1_ct*log_lam_safe )),0.5) ;
      log_lik_qpcr[r]  += normal_lpdf(qpcr_ct[r] | mu_ct, fmax(sigma_ct, 1e-6));
    }

    pp_qpcr_detect[r] = bernoulli_rng(p_det);
    if (pp_qpcr_detect[r] == 1) {
      real log_lam_safe = fmax(fmin(log_lam, 15.0), -10.0);
      real mu_ct        = alpha_ct - beta_ct * log_lam_safe;
      real sigma_ct     = pow(pow(sigma0_ct,2) + exp(2*(gamma0_ct + gamma1_ct*log_lam_safe )),0.5) ;
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

    vector[S+1] lam_edna_i;
    for (s in 1:(S+1)) lam_edna_i[s] = exp(fmin(log_lambda_edna[i, s], 15.0));

    // FIX 4: clamp before log
    real lam_sum     = fmax(sum(lam_edna_i), 1e-12);
    real log_lam_sum = log(lam_sum);

    vector[S+1] pi_i;
    if (lam_sum < 1e-12) {
      pi_i = rep_vector(1.0 / (S+1), S+1);
    } else {
      pi_i = lam_edna_i / lam_sum;
    }

    for (s in 1:(S+1)) {
      real log_lam_s  = fmax(fmin(log_lambda_edna[i, s], 15.0), -15.0);
      real lam_s      = exp(log_lam_s);

      real log_phi_s  = beta0_phi + log_lam_sum
                        + fmax(gamma0_phi - gamma1_phi * log_lam_s, 0.0);
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
