# =============================================================================
# 04_run_whale_edna_model_vS1.R
#
# Compile stan/whale_edna_hsgp_vS1.stan and fit it to the Stan data list
# built by scripts/03_format_stan_data_vS1.r.
#
# Inputs:  outputs/whale_edna_output_vS1/stan_data.rds
# Outputs: outputs/whale_edna_output_vS1/whale_edna_fit.rds
#          + CmdStan output files in outputs/whale_edna_output_vS1/
#
# =============================================================================

library(rstan)
# library(cmdstanr)
options(mc.cores = parallel::detectCores()) # rstan options.

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
# 1. Load pre-built Stan data
# -----------------------------------------------------------------------------
cat("=== Loading Stan data ===\n")
stan_data <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))

S       <- stan_data$S
M_total <- stan_data$M
cat(sprintf("  N=%d  S=%d  M=%d  N_qpcr=%d  N_mb=%d\n",
            stan_data$N, S, M_total, stan_data$N_qpcr_long, stan_data$N_mb_long))

### OLE CANNOT GET CMDSTANR to work
# # -----------------------------------------------------------------------------
# # 2. Compile
# # -----------------------------------------------------------------------------
# cat("=== Compiling Stan model ===\n")
# 
# mod <- cmdstan_model(
#   "stan/whale_edna_hsgp_vS1.stan",
#   cpp_options = list(stan_threads = TRUE),
#   force_recompile=TRUE
# )
# 
# init_fn <- function() {
#   list(
#     mu_sp      = c(runif(S,1.8,2),runif(1,4,4.1)),
#     gp_sigma   = runif(S,0.8,1.2), 
#     gp_l_raw   = matrix(runif(S*D1,0.5,0.6),
#                         nrow = S, ncol = D1, byrow = TRUE),
#     # gp_l       = runif(S*3,0.98,1.02)*matrix(c(50, 300, 150),
#     #                                          nrow = S, ncol = 3, byrow = TRUE),
#     z_beta     = matrix(0, nrow= S, ncol = M_total),
#     # log_RE_junk_raw = rnorm(N,0,1),
#     # sigma_junk = runif(1,0.2,0.4),
#     beta0_phi  = 0,
#     gamma0_phi = 2.0,
#     gamma1_phi = 1.0
#   )
# }
# # 
# # -----------------------------------------------------------------------------
# # 3. Sample
# # -----------------------------------------------------------------------------
# cat("=== Fitting model ===\n")
# fit <- mod$sample(
#   data              = stan_data,
#   chains            = N_CHAINS,
#   parallel_chains   = N_CHAINS,
#   threads_per_chain = 2,
#   iter_warmup       = N_WARMUP,
#   iter_sampling     = N_SAMPLE,
#   adapt_delta       = ADAPT_DELTA,
#   max_treedepth     = MAX_TREEDEPTH,
#   seed              = 42,
#   init              = init_fn,
#   output_dir        = OUTPUT_DIR,
#   show_messages     = TRUE
# )
# 

# THIS IS AN RSTAN VERSION BECAUSE OLE CAN'T GET
# THE cmdstanr package to my computer without updating R
# It is equivalent, perhaps not quite as good as the version Eiren 
# Originally had.
init_fn2 <- function(N_CHAIN){
  init_list <- list()
  for(i in 1:N_CHAIN){
    init_list[[i]] <- list(
      mu_sp      = c(runif(S,1.8,2),runif(1,4,4.1)),
      gp_sigma   = runif(S,0.8,1.2), 
      gp_l_raw   = matrix(runif(S*D1,0.2,0.3),
                          nrow = S, ncol = D1, byrow = TRUE),
      z_beta     = matrix(0, nrow= S, ncol = M_total),
      beta0_phi  = 0,
      gamma0_phi = 2.0,
      gamma1_phi = 1.0
    )
  }
  return(init_list)
}

stanMod = stan(file = "stan/whale_edna_hsgp_vS1.stan",
               data = stan_data, 
               verbose = FALSE, chains = N_CHAINS, thin = 1, 
               warmup = N_WARMUP, iter = N_WARMUP + N_SAMPLE, 
               init = init_fn2(N_CHAINS),
               control = list(max_treedepth=MAX_TREEDEPTH,
                              adapt_init_buffer = 75,
                              #                stepsize=0.01,
                              adapt_delta=ADAPT_DELTA),
               #                metric="diag_e"),
               #refresh = 100,
               sample_file="tmpX.csv",
               boost_lib = NULL
)

# - Load helper functions

source("scripts/helper_functions.R")


# Summarize fits. 
stanMod_summary <- summary(stanMod)$summary
pars <- extract(stanMod)

stanMod_main <- list()
## GP parameters.
gp_pars = c("lp__", "mu_sp","gp_sigma","gp_l_raw","gp_l")
stanMod_main$GP <- summary(stanMod,pars=gp_pars)$summary
traceplot(stanMod,pars=c("lp__",gp_pars))

## MB parameters
mb_pars <- c("beta0_phi", "gamma0_phi","gamma1_phi")
stanMod_main$MB <- summary(stanMod,pars=mb_pars)$summary
traceplot(stanMod,pars=mb_pars)

## QPCR parameters  
# qpcr_pars <- c("alpha_ct", "beta_ct","kappa","gamma0_ct","gamma1_ct","sigma0_ct")
# stanMod_main$MB <- summary(stanMod,pars=qpcr_pars)$summary
# traceplot(stanMod,pars=qpcr_pars)

# Summarize a few kinds of predictions and combine them with true
# ANIMALS
pred_obs_loc <- pred_summary(stanMod,
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
pred_all_loc <- pred_summary(stanMod,
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
pp    <- pred_summary(stanMod,pars=c("pp_pi_hat","pp_mb_reads"))


fits <- list(stanMod = stanMod,
             stan_data = stan_data,
             stanMod_main = stanMod_main,
             stanMod_summary = stanMod_summary,
             pars = pars, # posterior draws.
             
             pred_all_loc = pred_all_loc,
             pred_obs_loc = pred_obs_loc,
             pp = pp)

saveRDS(fits,file=file.path(OUTPUT_DIR, "whale_edna_vS1_fit.rds"))

