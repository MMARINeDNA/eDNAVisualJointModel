# =============================================================================
# 04_run_whale_edna_model_v4.2.R
#
# Compile stan/whale_edna_hsgp_v4.2.stan and fit it to the Stan data list
# built by scripts/03_format_stan_data_v4.2.r.
#
# Inputs:  outputs/whale_edna_output_v4.2/stan_data.rds
# Outputs: outputs/whale_edna_output_v4.2/whale_edna_fit.rds
#          + CmdStan output files in outputs/whale_edna_output_v4.2/
#
# v4.2 sampling configuration (per v3.2 lessons):
#   * iter_warmup = 1000, iter_sampling = 1000 (was 500 / 500)
#   * max_treedepth = 14 (was 12) - the M=3584 basis with stiff
#     posterior geometry needs longer trajectories.
#   * kappa, sigma_ct dropped from init_fn (now data, not parameters).
# =============================================================================

#library(cmdstanr)
library(rstan)
options(mc.cores = parallel::detectCores()) # rstan options.
set.seed(42)

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
N_CHAINS      <- 3
N_WARMUP      <- 200
N_SAMPLE      <- 100
ADAPT_DELTA   <- 0.80
MAX_TREEDEPTH <- 12

OUTPUT_DIR <- "outputs/whale_edna_output_v4.2"
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

 # THIS IS AN RSTAN VERSION BECAUSE OLE CAN"T PULL
 # THE cmdstanr package to my computer without updating R
 # It is equivalent, perhaps not quite as good as the version Eiren 
 # Originally had.
init_fn2 <- function(N_CHAIN){
  init_list <- list()
  for(i in 1:N_CHAIN){
    init_list[[i]] <- list(
      mu_sp      = runif(S+1,1.8,2),
      gp_sigma   = runif(S,0.8,1.2),
      gp_l       = runif(S*3,0.98,1.02)*matrix(c(50, 300, 150),
                          nrow = S, ncol = 3, byrow = TRUE),
      z_beta     = matrix(0, nrow = S, ncol = M_total),
      log_RE_junk_raw = runif(N,-0.01,0.01),
      sigma_junk = runif(1,0.5,1.0),
      beta0_phi  = 0,
      gamma0_phi = 2.0,
      gamma1_phi = 1.0
    )
  }
  return(init_list)
}

stanMod = stan(file = "stan/whale_edna_hsgp_v4.2.stan",
               data = stan_data, 
               verbose = FALSE, chains = N_CHAINS, thin = 1, 
               warmup = N_WARMUP, iter = N_WARMUP + N_SAMPLE, 
               init = init_fn2(N_CHAINS),
               control = list(max_treedepth=MAX_TREEDEPTH,
                              adapt_init_buffer = 75,
               #                stepsize=0.01,
                               adapt_delta=0.8),
               #                metric="diag_e"),
               #refresh = 100,
               sample_file="tmpY.csv",
               boost_lib = NULL
)






#### EIRENS ORIGINAL STAN CALL.


init_fn <- function(nchain) {
  list(
    mu_sp      = rep(2.0, S+1),
    gp_sigma   = rep(1.0, S),
    gp_l       = matrix(c(50, 300, 150),
                        nrow = S, ncol = 3, byrow = TRUE),
    z_beta     = matrix(0.0, nrow = S, ncol = M_total),
    beta0_phi  = 2.0,
    gamma0_phi = 5.0,
    gamma1_phi = 1.0
  )
}

mod <- cmdstan_model(
  "stan/whale_edna_hsgp_v4.2.stan",
  cpp_options = list(stan_threads = TRUE)
)

init_fn <- function() {
  list(
    mu_sp      = rep(2.0, S+1),
    gp_sigma   = rep(1.0, S),
    gp_l       = matrix(c(50, 300, 150),
                        nrow = S, ncol = 3, byrow = TRUE),
    z_beta     = matrix(0.0, nrow = S, ncol = M_total),
    beta0_phi  = 2.0,
    gamma0_phi = 2.0,
    gamma1_phi = 1.0
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
