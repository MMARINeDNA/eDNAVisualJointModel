// =============================================================================
// qpcr_standards.stan
//
// Standard curves for the hake qPCR, one per plate, fitted to the gBlocks
// standards. After the qPCR vignette's qPCR_standards_multi-plate.stan
// (Shelton et al.), with its priors and its bound on the slope.
//
// For standard i on plate p with K_i copies per uL and a 2 uL reaction
// (lambda_i = 2 K_i expected copies in the well):
//   detection   G_i ~ Bernoulli(1 - exp(-kappa * lambda_i))
//   Ct | G_i=1  C_i ~ Normal(beta0_p + beta1_p * log(lambda_i), sigma_C(lambda_i))
//   sigma_C(lambda) = sqrt(sigma0^2 + exp(2 (gamma0 + gamma1 log lambda)))
// beta0_p and beta1_p are per plate; sigma0, gamma0, gamma1 and kappa are
// shared. beta1_p <= -1/log(2): copies can at most double each cycle.
//
// In the eDNA model's qPCR hurdle (HSGP4eDNA edna_functions.stan):
//   alpha_ct = beta0_p, beta_ct = -beta1_p, kappa = kappa,
//   gamma0_ct = gamma0, gamma1_ct = gamma1, sigma0_ct = sigma0,
// with lambda = expected copies in the reaction.
// =============================================================================

data {
  int<lower=1> N;                          // standards (wells)
  int<lower=1> N_plate;
  array[N] int<lower=1, upper=N_plate> plate;
  array[N] int<lower=0, upper=1> detect;
  vector[N] ct;                            // Ct; ignored where detect = 0
  vector<lower=0>[N] K;                    // copies per uL
  real<lower=0> vol_reaction;              // uL of template per reaction (2)

  vector[2] beta0_prior;                   // normal(mean, sd)
  vector[2] beta1_prior;
  real beta1_max;                          // -1 / log(2)
  vector[2] gamma0_prior;
  vector[2] gamma1_prior;
  vector[2] logit_kappa_prior;
}

transformed data {
  vector[N] log_lam = log(vol_reaction * K);
  int N_pos = sum(detect);
  array[N_pos] int pos;
  {
    int j = 0;
    for (i in 1:N) if (detect[i] == 1) { j += 1; pos[j] = i; }
  }
}

parameters {
  vector[N_plate] beta0;
  vector<upper=beta1_max>[N_plate] beta1;
  real<lower=0> sigma0;
  real gamma0;
  real gamma1;
  real logit_kappa;
}

transformed parameters {
  real kappa = inv_logit(logit_kappa);
}

model {
  beta0 ~ normal(beta0_prior[1], beta0_prior[2]);
  beta1 ~ normal(beta1_prior[1], beta1_prior[2]);
  sigma0 ~ normal(0, 1);
  gamma0 ~ normal(gamma0_prior[1], gamma0_prior[2]);
  gamma1 ~ normal(gamma1_prior[1], gamma1_prior[2]);
  logit_kappa ~ normal(logit_kappa_prior[1], logit_kappa_prior[2]);

  detect ~ bernoulli(1 - exp(-kappa * exp(log_lam)));
  {
    vector[N_pos] lp = log_lam[pos];
    ct[pos] ~ normal(beta0[plate[pos]] + beta1[plate[pos]] .* lp,
                     sqrt(square(sigma0) + exp(2 * (gamma0 + gamma1 * lp))));
  }
}

generated quantities {
  // In the eDNA model's terms
  vector[N_plate] alpha_ct = beta0;
  vector[N_plate] beta_ct  = -beta1;
  vector[N_plate] efficiency = exp(-1 ./ beta1) - 1;   // fraction more copies per cycle (1 = doubling)
  // Posterior predictive Ct for every standard (for detected ones, the check)
  vector[N] ct_rep;
  for (i in 1:N)
    ct_rep[i] = normal_rng(beta0[plate[i]] + beta1[plate[i]] * log_lam[i],
                           sqrt(square(sigma0) + exp(2 * (gamma0 + gamma1 * log_lam[i]))));
}
