# =============================================================================
# 04_run_whale_edna_model_vS1.R
#
# Compile stan/whale_edna_hsgp_vS1.stan and fit it to the Stan data list
# built by scripts/simulation/03_format_stan_data_vS1.r.
#
# Inputs:  outputs/whale_edna_output_vS1/stan_data.rds
#          outputs/whale_edna_output_vS1/pred_points.rds
#          outputs/whale_edna_output_vS1/whale_edna_sim_vS1.rds
# Outputs: outputs/whale_edna_output_vS1/whale_edna_vS1_fit.rds
#          + CmdStan CSV output files in outputs/whale_edna_output_vS1/
#
# Runner: cmdstanr is the canonical sampler (compile + sample) and all
# post-processing uses the cmdstanr / posterior API (fit$draws(),
# fit$summary(), fit$sampler_diagnostics()). An earlier interim version
# sampled with rstan::stan() because cmdstanr was not installable on one
# machine; that path has been removed now that cmdstanr is the standard.
#
# The saved fit bundle stores only plain serialisable objects (posterior
# draws arrays, summary tibbles, prediction data frames) - not the live
# CmdStanMCMC object - so 05 can load it standalone.
#
# This script is self-contained: it loads everything it needs from disk and
# can be run on its own after 03 has written the outputs above.
# =============================================================================

library(cmdstanr)   # canonical sampler (compile + sample)
library(posterior)  # draws summaries
library(tidyverse)

options(mc.cores = parallel::detectCores())
set.seed(42)

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
N_CHAINS      <- 3
N_WARMUP      <- 400
N_SAMPLE      <- 600
ADAPT_DELTA   <- 0.70
MAX_TREEDEPTH <- 11

OUTPUT_DIR <- "outputs/whale_edna_output_vS1"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load pre-built Stan data, prediction grid, and simulation truth
# -----------------------------------------------------------------------------
cat("=== Loading Stan data / prediction grid / simulation ===\n")
stan_data   <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))
pred_points <- readRDS(file.path(OUTPUT_DIR, "pred_points.rds"))
sim         <- readRDS(file.path(OUTPUT_DIR, "whale_edna_sim_vS1.rds"))

S       <- stan_data$S
M_total <- stan_data$M
D1      <- stan_data$D1
cat(sprintf("  N=%d  S=%d  M=%d  D1=%d  N_qpcr=%d  N_mb=%d  N_pred=%d\n",
            stan_data$N, S, M_total, D1,
            stan_data$N_qpcr_long, stan_data$N_mb_long, stan_data$N_pred))

# -----------------------------------------------------------------------------
# 2. Compile (cmdstanr)
# -----------------------------------------------------------------------------
cat("=== Compiling Stan model (cmdstanr) ===\n")
mod <- cmdstan_model("stan/whale_edna_hsgp_vS1.stan")

# One init list per chain (cmdstanr calls this once per chain).
init_fn <- function() {
  list(
    mu_sp      = c(runif(S, 1.8, 2.0), runif(1, 4.0, 4.1)),
    gp_sigma   = runif(S, 0.8, 1.2),
    gp_l_raw   = matrix(runif(S * D1, 0.2, 0.3),
                        nrow = S, ncol = D1, byrow = TRUE),
    z_beta     = matrix(0, nrow = S, ncol = M_total),
    beta0_phi  = 0,
    gamma0_phi = 2.0,
    gamma1_phi = 1.0
  )
}

# -----------------------------------------------------------------------------
# 3. Sample (cmdstanr)
# -----------------------------------------------------------------------------
cat("=== Fitting model ===\n")
fit_cmd <- mod$sample(
  data            = stan_data,
  chains          = N_CHAINS,
  parallel_chains = N_CHAINS,
  iter_warmup     = N_WARMUP,
  iter_sampling   = N_SAMPLE,
  adapt_delta     = ADAPT_DELTA,
  max_treedepth   = MAX_TREEDEPTH,
  seed            = 42,
  init            = init_fn,
  output_dir      = OUTPUT_DIR,
  show_messages   = TRUE
)

# -----------------------------------------------------------------------------
# 4. Summaries and predictions (cmdstanr / posterior API)
# -----------------------------------------------------------------------------
source("scripts/helper_functions.R")

# Model-parameter draws only (NOT the large generated-quantity arrays) — this
# is what 05 needs for diagnostics, trace plots, and parameter recovery, and
# it keeps the saved bundle small.
param_pars <- c("lp__", "mu_sp", "gp_sigma", "gp_l_raw", "gp_l",
                "beta0_phi", "gamma0_phi", "gamma1_phi")
param_draws  <- fit_cmd$draws(variables = param_pars)
sampler_diag <- fit_cmd$sampler_diagnostics()
diag_summary <- fit_cmd$diagnostic_summary()

# Per-group summary tibbles (variable, mean, sd, rhat, ess_bulk, ess_tail, ...).
gp_pars <- c("lp__", "mu_sp", "gp_sigma", "gp_l_raw", "gp_l")
mb_pars <- c("beta0_phi", "gamma0_phi", "gamma1_phi")
stanMod_main <- list(
  GP = fit_cmd$summary(variables = gp_pars),
  MB = fit_cmd$summary(variables = mb_pars)
)
stanMod_summary <- dplyr::bind_rows(stanMod_main$GP, stanMod_main$MB)

## QPCR parameters are fixed as data in vS1 (not sampled).

# Summarize a few kinds of predictions and combine them with true
# ANIMALS
pred_obs_loc <- pred_summary(fit_cmd,
                      pars=c("log_lambda","log_lambda_edna"))
pred_obs_loc$samples <- sim$design$samples

pred_obs_loc$log_lambda$Long <- left_join(pred_obs_loc$log_lambda$Long,
                                    pred_obs_loc$samples,by = join_by(idx == station))
A           <- sim$truth$lambda_true_si %>% as.data.frame()
colnames(A) <- paste0("sp_",1:ncol(A))
A$idx       <- 1:nrow(A)
B<- pivot_longer(A,cols= -idx,names_to = "species",values_to = "true_lambda")  %>%
      mutate(true_log_lambda = log(true_lambda))
pred_obs_loc$log_lambda$Long <- left_join(pred_obs_loc$log_lambda$Long, B) %>%
                                  mutate(resid = Mean - true_log_lambda)

pred_obs_loc$log_lambda_edna$Long <- left_join(pred_obs_loc$log_lambda_edna$Long,
                                          pred_obs_loc$samples,by = join_by(idx == station))
A           <- sim$truth$C_obs_si %>% as.data.frame()
colnames(A) <- paste0("sp_",1:ncol(A))
A$idx       <- 1:nrow(A)
B <- pivot_longer(A,cols= -idx,names_to = "species",values_to = "true_lambda_edna")  %>%
        mutate(true_log_lambda_edna = log(true_lambda_edna))
pred_obs_loc$log_lambda_edna$Long <- left_join(pred_obs_loc$log_lambda_edna$Long, B) %>%
                                    mutate(resid = Mean - true_log_lambda_edna)

## EDNA COPIES
pred_all_loc <- pred_summary(fit_cmd,
                      pars=c("log_lambda_pred","log_lambda_edna_pred"))
pred_point_km <- pred_points
pred_point_km$idx <- 1:nrow(pred_point_km)
colnames(pred_point_km) <- c("X","Y","idx")
pred_all_loc$samples <- pred_point_km
pred_all_loc$pred_location <- cbind(idx=1:nrow(pred_points),pred_points)
pred_all_loc$log_lambda_pred$Long <- left_join(pred_all_loc$log_lambda_pred$Long,
                                          pred_all_loc$pred_location)
pred_all_loc$log_lambda_edna_pred$Long <- left_join(pred_all_loc$log_lambda_edna_pred$Long,
                                               pred_all_loc$pred_location)

## Posterior Predictions
pp    <- pred_summary(fit_cmd,pars=c("pp_pi_hat","pp_mb_reads"))


fits <- list(draws           = param_draws,    # posterior draws (model params)
             sampler_diag    = sampler_diag,    # NUTS sampler diagnostics
             diag_summary    = diag_summary,    # divergences / treedepth / ebfmi
             stan_data       = stan_data,
             stanMod_main    = stanMod_main,    # per-group summary tibbles
             stanMod_summary = stanMod_summary,

             pred_all_loc = pred_all_loc,
             pred_obs_loc = pred_obs_loc,
             pp = pp)

saveRDS(fits,file=file.path(OUTPUT_DIR, "whale_edna_vS1_fit.rds"))
cat(sprintf("Fit saved to %s\n", file.path(OUTPUT_DIR, "whale_edna_vS1_fit.rds")))
