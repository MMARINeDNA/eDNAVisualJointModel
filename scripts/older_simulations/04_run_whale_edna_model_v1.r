# =============================================================================
# 04_run_whale_edna_model_v1.R
#
# Compile stan/whale_edna_hsgp_v1.stan and fit it to the Stan data list built
# by scripts/03_format_stan_data_v1.r.
#
# Inputs:  outputs/whale_edna_output_v1/stan_data.rds
# Outputs: outputs/whale_edna_output_v1/whale_edna_fit.rds
# =============================================================================

library(cmdstanr)

set.seed(42)

N_CHAINS      <- 4
N_WARMUP      <- 250
N_SAMPLE      <- 250
ADAPT_DELTA   <- 0.90
MAX_TREEDEPTH <- 12

OUTPUT_DIR <- "outputs/whale_edna_output_v1"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Loading Stan data ===\n")
stan_data <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))

S       <- stan_data$S
M_total <- stan_data$M
cat(sprintf("  N=%d  S=%d  R=%d  K=%d  M=%d\n",
            stan_data$N, S, stan_data$R, stan_data$K, M_total))

cat("=== Compiling Stan model ===\n")
mod <- cmdstan_model("stan/whale_edna_hsgp_v1.stan")

init_fn <- function() {
  list(
    gp_sigma = rep(0.5, S),
    gp_l     = matrix(c(50, 150, 300), nrow = S, ncol = 3, byrow = TRUE),
    z_beta   = matrix(0, nrow = S, ncol = M_total),
    mu_sp    = rep(2.0, S),
    kappa    = 0.8,
    alpha_ct = 38.0,
    beta_ct  = 3.32,
    sigma_ct = 0.5,
    p_zi     = 0.3,
    phi_bb   = rep(10.0, S)
  )
}

cat("=== Fitting model ===\n")
fit <- mod$sample(
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

fit$save_object(file.path(OUTPUT_DIR, "whale_edna_fit.rds"))
cat(sprintf("Fit saved to %s\n", file.path(OUTPUT_DIR, "whale_edna_fit.rds")))
