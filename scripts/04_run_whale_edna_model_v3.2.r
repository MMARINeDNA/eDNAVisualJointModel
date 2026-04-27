# =============================================================================
# 04_run_whale_edna_model_v3.2.R
#
# Compile stan/whale_edna_hsgp_v3.2.stan and fit it to the Stan data list
# built by scripts/03_format_stan_data_v3.2.r.
#
# Inputs:  outputs/whale_edna_output_v3.2/stan_data.rds
# Outputs: outputs/whale_edna_output_v3.2/whale_edna_fit.rds
#          + CmdStan output files in outputs/whale_edna_output_v3.2/
#
# v3.2: kappa and sigma_ct are fixed in data; metabarcoding
# overdispersion params (beta0_phi / gamma0_phi / gamma1_phi) are gone.
# Only mu_sp, gp_sigma, gp_l, z_beta are sampled.
# =============================================================================

library(cmdstanr)

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
N_CHAINS      <- 4
N_WARMUP      <- 1000   # bumped from 500: with M=3584 and a near-degenerate
                        # gp_sigma <-> z_beta ridge in the posterior, NUTS
                        # needs more warmup to find the correct mode and
                        # adapt away from the low-sigma trap.
N_SAMPLE      <- 500
ADAPT_DELTA   <- 0.90
MAX_TREEDEPTH <- 12

OUTPUT_DIR <- "outputs/whale_edna_output_v3.2"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load pre-built Stan data
# -----------------------------------------------------------------------------
cat("=== Loading Stan data ===\n")
stan_data <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))

S       <- stan_data$S
M_total <- stan_data$M
cat(sprintf("  N=%d  S=%d  M=%d  N_qpcr=%d\n",
            stan_data$N, S, M_total, stan_data$N_qpcr_long))

# -----------------------------------------------------------------------------
# 2. Compile
# -----------------------------------------------------------------------------
cat("=== Compiling Stan model ===\n")

mod <- cmdstan_model(
  "stan/whale_edna_hsgp_v3.2.stan",
  cpp_options = list(stan_threads = TRUE)
)

init_fn <- function() {
  list(
    mu_sp      = rep(2.0, S),
    gp_sigma   = rep(0.8, S),
    gp_l       = matrix(c(50, 300, 150),
                        nrow = S, ncol = 3, byrow = TRUE),
    z_beta     = matrix(0.0, nrow = S, ncol = M_total)
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
