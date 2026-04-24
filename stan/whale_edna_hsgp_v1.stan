// =============================================================================
// whale_edna_hsgp.stan
//
// 3-D anisotropic HSGP model for three-species whale/fish eDNA.
// Oregon / Washington coast, UTM Zone 10N.
//
// Species:
//   s=1  Merluccius productus        (hake)   — qPCR + metabarcoding
//   s=2  Megaptera novaeangliae      (humpback) — metabarcoding only
//   s=3  Lagenorhynchus obliquidens  (PWSD)     — metabarcoding only
//
// Latent field:
//   log(lambda_si) = mu_s
//                  + f_s(X_i, Y_i, Z_bathy_i)    [HSGP over bathymetry]
//                  + log_zsample_effect[i, s]     [fixed water column offset]
//
//   f_s ~ HSGP(sigma_s, [lx_s, ly_s, lz_s])
//
// Coordinates passed in are normalised (X, Y, Z_bathy) — Z_bathy is the
// bottom depth at each station, NOT the sample collection depth.
// The water column effect (log_zsample_effect) is pre-computed in R and
// passed as fixed data.
//
// Observation model 1 — qPCR hurdle (hake only, R=3 replicates):
//   p_det_si      = 1 - exp(-kappa * lambda_edna_si * vol_i)
//   any_detect_si ~ Bernoulli(1 - (1 - p_det_si)^R)
//   mean_ct_si | detect ~ Normal(alpha_ct - beta_ct*log(lambda_edna_si*vol_i),
//                                sigma_ct / sqrt(n_detect_si))
//
// Observation model 2 — metabarcoding Beta-Binomial (all species):
//   MARVER3 marker; K=3 replicates per sample, each modelled independently.
//   pi_si         = lambda_edna_si / sum_s(lambda_edna_si)
//   mb_reads_sik  ~ BetaBinomial(mb_total_ik, mu=pi_si, phi=phi_bb_s)
//
// HSGP basis: Riutort-Mayol et al. (2023), Statistics & Computing.
// =============================================================================

functions {

  // -----------------------------------------------------------------------
  // 1-D Laplacian eigenfunction on [-L, L]
  // -----------------------------------------------------------------------
  real phi1d(real L, int m, real x) {
    return (1.0 / sqrt(L)) * sin(m * pi() * (x + L) / (2.0 * L));
  }

  // -----------------------------------------------------------------------
  // Spectral density of 1-D squared-exponential kernel
  // -----------------------------------------------------------------------
  real spd_se_1d(real alpha, real rho, real w) {
    return square(alpha) * sqrt(2.0 * pi()) * rho
             * exp(-0.5 * square(rho * w));
  }

  // -----------------------------------------------------------------------
  // Build N x M HSGP design matrix (3-D)
  // -----------------------------------------------------------------------
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

  // -----------------------------------------------------------------------
  // Spectral weights (length M) for anisotropic 3-D SE kernel
  // -----------------------------------------------------------------------
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

  // -----------------------------------------------------------------------
  // Beta-Binomial log PMF (mean / precision parameterisation)
  // Guards against degenerate proportions pushing alpha or beta to zero.
  // -----------------------------------------------------------------------
  real beta_binomial2_lpmf(int y, int n, real mu, real phi) {
    real mu_safe  = fmin(fmax(mu,  1e-6), 1.0 - 1e-6);
    real phi_safe = fmax(phi, 1e-6);
    return beta_binomial_lpmf(y | n, mu_safe * phi_safe,
                                     (1.0 - mu_safe) * phi_safe);
  }

}

// =============================================================================
data {

  int<lower=1> N;    // total samples
  int<lower=1> S;    // species (3)
  int<lower=1> R;    // qPCR replicates (hake only; = 3)
  int<lower=1> K;    // metabarcoding replicates per sample (= 3)
  int<lower=1> M;    // HSGP basis dimension = m1*m2*m3

  // Spatial coordinates normalised to [-1, 1]
  // col 1 = X (km), col 2 = Y (km), col 3 = Z_bathy (m)
  matrix[N, 3] coords;

  // Half-ranges used for normalisation (original units)
  vector<lower=0>[3] coord_scale;

  // HSGP configuration
  vector<lower=0>[3]  L_hsgp;
  array[3] int<lower=1> m_hsgp;

  // Sample volumes (litres) — must be > 0
  vector<lower=1e-6>[N] vol;

  // Water column eDNA effect (fixed, pre-computed in R)
  // log_zsample_effect[i, s] = log multiplier on lambda at sample depth
  matrix[N, S] log_zsample_effect;

  // --- qPCR data (hake only, s=1) ---
  array[N] int<lower=0, upper=R> n_detect;
  vector[N]                      mean_ct_obs;   // 0 if n_detect == 0
  array[N] int<lower=0, upper=1> any_detect;

  // --- Metabarcoding data (all species, K replicates) ---
  array[N, S, K] int<lower=0> mb_reads;
  array[N, K]    int<lower=0> mb_total;

}

// =============================================================================
transformed data {

  // Design matrix shared across species
  matrix[N, M] PHI = hsgp_phi(coords, L_hsgp, m_hsgp);

  // Log-volume: guarded against zero (vol is constrained > 1e-6 in data)
  vector[N] log_vol;
  for (i in 1:N) {
    log_vol[i] = log(fmax(vol[i], 1e-6));
  }

}

// =============================================================================
parameters {

  // --- GP hyperparameters (per species) ---
  vector<lower=0>[S]    gp_sigma;
  matrix<lower=0>[S, 3] gp_l;       // [species × dim]: lx(km), ly(km), lz(m)

  // --- Non-centred GP basis coefficients ---
  matrix[S, M] z_beta;

  // --- Log-scale species intercepts ---
  vector[S] mu_sp;

  // --- qPCR observation parameters (hake only) ---
  real<lower=0>          kappa;      // detection sensitivity
  // Bounded to guarantee mu_ct is always finite:
  // mu_ct = alpha_ct - beta_ct * log_lam_clamped
  // log_lam_clamped in [-10, 15], so with these bounds
  // mu_ct in [20 - 10*15, 50 - 0*(-10)] = [-100, 50]: always finite
  real<lower=20, upper=50> alpha_ct;
  real<lower=0,  upper=10> beta_ct;
  real<lower=0>            sigma_ct;

  // --- Zero-inflation probability (hake only) ---
  real<lower=0, upper=1> p_zi;

  // --- Metabarcoding overdispersion (per species) ---
  vector<lower=0.1>[S] phi_bb;

}

// =============================================================================
transformed parameters {

  // log(lambda_si): log true animal density — function of (X, Y, Z_bathy)
  matrix[N, S] log_lambda;

  // log(lambda_edna_si): log eDNA concentration at sample collection depth
  matrix[N, S] log_lambda_edna;

  for (s in 1:S) {
    vector[M] wt = hsgp_weights(
      gp_sigma[s],
      to_vector(gp_l[s, ]),
      L_hsgp,
      coord_scale,
      m_hsgp
    );
    vector[N] f_s = PHI * (wt .* to_vector(z_beta[s, ]));

    for (i in 1:N) {
      log_lambda[i, s]      = mu_sp[s] + f_s[i];
      log_lambda_edna[i, s] = log_lambda[i, s] + log_zsample_effect[i, s];
    }
  }

}

// =============================================================================
model {

  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------

  gp_sigma ~ normal(0, 1.5);

  for (s in 1:S) {
    gp_l[s, 1] ~ normal(50,  40);    // lx (km)
    gp_l[s, 2] ~ normal(150, 80);    // ly (km)
    gp_l[s, 3] ~ normal(300, 150);   // lz_bathy (m)
  }

  to_vector(z_beta) ~ std_normal();

  mu_sp ~ normal(2.0, 1.5);

  // qPCR — priors consistent with parameter bounds
  kappa    ~ normal(1.0, 0.5);
  alpha_ct ~ normal(38.0, 3.0);
  beta_ct  ~ normal(3.32, 0.5);
  sigma_ct ~ normal(0.5,  0.3);
  p_zi     ~ beta(1.5, 4.0);

  phi_bb ~ normal(10.0, 5.0);

  // ------------------------------------------------------------------
  // Likelihood
  // ------------------------------------------------------------------
  for (i in 1:N) {

    // eDNA lambda vector at sample depth
    vector[S] lam_edna_i;
    for (s in 1:S) lam_edna_i[s] = exp(log_lambda_edna[i, s]);

    // Compositional proportions — clamped and renormalised
    vector[S] pi_i = lam_edna_i / sum(lam_edna_i);
    for (s in 1:S) pi_i[s] = fmin(fmax(pi_i[s], 1e-6), 1.0 - 1e-6);
    pi_i = pi_i / sum(pi_i);

    // ---- qPCR hurdle (hake, s=1) ------------------------------------
    {
      // Single clamped quantity used for ALL qPCR arithmetic
      real log_lam_clamped = fmax(fmin(
        log_lambda_edna[i, 1] + log_vol[i], 15.0), -10.0);
      real lam_edna_vol = exp(log_lam_clamped);

      real p_det_rep = 1.0 - exp(-kappa * lam_edna_vol);
      real p_eff     = (1.0 - p_zi) * p_det_rep;
      real p_any     = 1.0 - pow(1.0 - p_eff, R);

      target += bernoulli_lpmf(any_detect[i] | p_any);

      if (any_detect[i] == 1 && n_detect[i] > 0) {
        real mu_ct = alpha_ct - beta_ct * log_lam_clamped;
        real sd_ct = fmax(sigma_ct / sqrt(n_detect[i] + 0.0), 1e-6);
        target += normal_lpdf(mean_ct_obs[i] | mu_ct, sd_ct);
      }
    }

    // ---- Metabarcoding Beta-Binomial (all species, K replicates) ----
    for (k in 1:K) {
      if (mb_total[i, k] > 0) {
        for (s in 1:S) {
          target += beta_binomial2_lpmf(
            mb_reads[i, s, k] | mb_total[i, k],
            pi_i[s], phi_bb[s]
          );
        }
      }
    }

  }
}

// =============================================================================
generated quantities {

  // Log-likelihoods for LOO-CV
  vector[N]    log_lik_qpcr;
  matrix[N, S] log_lik_mb;

  // Posterior predictive
  array[N] int   pp_n_detect;
  vector[N]      pp_mean_ct;
  array[N, S, K] int pp_mb_reads;

  // Posterior mean densities
  matrix[N, S] lambda_hat;       // true animal density (no water column effect)
  matrix[N, S] lambda_edna_hat;  // eDNA at sample depth

  for (i in 1:N) {

    // Build eDNA lambda and proportions
    vector[S] lam_edna_i;
    for (s in 1:S) {
      lam_edna_i[s]         = exp(log_lambda_edna[i, s]);
      lambda_hat[i, s]      = exp(log_lambda[i, s]);
      lambda_edna_hat[i, s] = lam_edna_i[s];
    }
    vector[S] pi_i = lam_edna_i / sum(lam_edna_i);
    for (s in 1:S) pi_i[s] = fmin(fmax(pi_i[s], 1e-6), 1.0 - 1e-6);
    pi_i = pi_i / sum(pi_i);

    // ---- qPCR (hake) ------------------------------------------------
    {
      real log_lam_clamped = fmax(fmin(
        log_lambda_edna[i, 1] + log_vol[i], 15.0), -10.0);
      real lam_edna_vol = exp(log_lam_clamped);
      real p_det_rep    = 1.0 - exp(-kappa * lam_edna_vol);
      real p_eff        = (1.0 - p_zi) * p_det_rep;
      real p_any        = 1.0 - pow(1.0 - p_eff, R);

      // Log-likelihood
      log_lik_qpcr[i] = bernoulli_lpmf(any_detect[i] | p_any);
      if (any_detect[i] == 1 && n_detect[i] > 0) {
        real mu_ct = alpha_ct - beta_ct * log_lam_clamped;
        real sd_ct = fmax(sigma_ct / sqrt(n_detect[i] + 0.0), 1e-6);
        log_lik_qpcr[i] += normal_lpdf(mean_ct_obs[i] | mu_ct, sd_ct);
      }

      // Posterior predictive
      {
        int  any_pp = bernoulli_rng(p_any);
        int  nd     = 0;
        real ct_pp  = 0.0;
        if (any_pp == 1) {
          for (r in 1:R) nd += bernoulli_rng(p_eff);
          nd = max(nd, 1);
          real mu_ct = alpha_ct - beta_ct * log_lam_clamped;
          ct_pp = normal_rng(mu_ct,
                             fmax(sigma_ct / sqrt(nd + 0.0), 1e-6));
        }
        pp_n_detect[i] = nd;
        pp_mean_ct[i]  = ct_pp;
      }
    }

    // ---- Metabarcoding (all species, K replicates) ------------------
    for (s in 1:S) log_lik_mb[i, s] = 0.0;

    for (k in 1:K) {
      for (s in 1:S) {
        if (mb_total[i, k] > 0) {
          log_lik_mb[i, s] += beta_binomial2_lpmf(
            mb_reads[i, s, k] | mb_total[i, k],
            pi_i[s], phi_bb[s]
          );
          real a_bb = pi_i[s] * phi_bb[s];
          real b_bb = (1.0 - pi_i[s]) * phi_bb[s];
          real p_bb = beta_rng(a_bb, b_bb);
          pp_mb_reads[i, s, k] = binomial_rng(mb_total[i, k], p_bb);
        } else {
          log_lik_mb[i, s]    += 0.0;
          pp_mb_reads[i, s, k] = 0;
        }
      }
    }

  }
}
