# =============================================================================
# 04_run_whale_edna_model_v4.R
#
# Compile stan/whale_edna_hsgp_v4.stan and fit it to the Stan data list built
# by scripts/03_format_stan_data_v4.r.
#
# Inputs:  outputs/whale_edna_output_v4/stan_data.rds
# Outputs: outputs/whale_edna_output_v4/whale_edna_fit.rds
#          + CmdStan output files in outputs/whale_edna_output_v4/
# =============================================================================

library(cmdstanr)

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
N_CHAINS      <- 4
N_WARMUP      <- 500
N_SAMPLE      <- 500
ADAPT_DELTA   <- 0.90
MAX_TREEDEPTH <- 12

OUTPUT_DIR <- "outputs/whale_edna_output_v4"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load pre-built Stan data
# -----------------------------------------------------------------------------
cat("=== Loading Stan data ===\n")
stan_data <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))

S       <- stan_data$S
M_total <- stan_data$M
cat(sprintf("  N=%d  S=%d  M=%d  N_qpcr=%d  N_mb=%d\n",
            stan_data$N, S, M_total, stan_data$N_qpcr_long, stan_data$N_mb_long))

# -----------------------------------------------------------------------------
# 2. Compile
# -----------------------------------------------------------------------------
cat("=== Compiling Stan model ===\n")

mod <- cmdstan_model(
  "stan/whale_edna_hsgp_v4.stan",
  cpp_options = list(stan_threads = TRUE)
)

init_fn <- function() {
  list(
    mu_sp      = rep(2.0, S),
    gp_sigma   = rep(0.8, S),
    gp_l       = matrix(c(50, 150, 300),
                        nrow = S, ncol = 3, byrow = TRUE),
    z_beta     = matrix(0.0, nrow = S, ncol = M_total),
    kappa      = 0.85,
    sigma_ct   = 0.50,
    beta0_phi  = rep(2.0, S),
    gamma0_phi = rep(5.0, S),
    gamma1_phi = rep(1.0, S)
  )
}

# -----------------------------------------------------------------------------
# 3. Sample
# -----------------------------------------------------------------------------
cat("=== Fitting model ===\n")
fit <- mod$sample(
  data              = stan_data,
  chains            = N_CHAINS,
  parallel_chains   = N_CHAINS,
  threads_per_chain = 2,
  iter_warmup       = N_WARMUP,
  iter_sampling     = N_SAMPLE,
  adapt_delta       = ADAPT_DELTA,
  max_treedepth     = MAX_TREEDEPTH,
  seed              = 42,
  init              = init_fn,
  output_dir        = OUTPUT_DIR,
  show_messages     = TRUE
)

fit$save_object(file.path(OUTPUT_DIR, "whale_edna_fit.rds"))
cat(sprintf("Fit saved to %s\n", file.path(OUTPUT_DIR, "whale_edna_fit.rds")))
