# =============================================================================
# 06_basis_sweep_vS1.R
#
# Systematic HSGP basis-count sensitivity sweep for the vS1 scenario.
#
# For each (basis config x simulation seed): simulate a vS1 dataset, fit the
# canonical model at that basis count, and record sampling health, GP
# hyperparameter recovery (lx, ly, gp_sigma, mu), and latent-field recovery
# (posterior-mean log_lambda vs truth at stations). Sweeping the basis from
# under- to over-resolved brackets both the oversmoothing bias (too few
# basis) and the overfitting signature (gp_l collapsing toward 0).
#
# Faithful representation rule (Riutort-Mayol): with boundary c = 1.5 and
# normalised length scales (lx, ly) = (0.20, 0.47), m >= (13, 6). Configs
# below bracket that. Configs H/I hold M fixed but swap the shape to show the
# short cross-shore scale is the binding constraint.
#
# The prediction grid is deliberately tiny here (predictions aren't used) so
# fits stay fast and CmdStan CSVs stay small. One compiled model serves all
# configs (basis count is data, not model).
#
# Results (incremental, resumable) -> outputs/whale_edna_output_vS1/basis_sweep/
#   sweep_results.rds : list(diag, recovery, field) of tidy long tibbles
#
# Env overrides (for quick smoke runs), all optional:
#   SWEEP_CONFIGS="B,E"  SWEEP_SEEDS="101"  SWEEP_CHAINS=2
#   SWEEP_WARMUP=40  SWEEP_SAMPLE=40  SWEEP_PREDN=4
# =============================================================================

library(cmdstanr)
library(posterior)
suppressMessages(library(tidyverse))
source("scripts/simulation/vS1_functions.R")

# -----------------------------------------------------------------------------
# 0. Sweep configuration
# -----------------------------------------------------------------------------
CONFIGS <- tribble(
  ~label, ~mx, ~my,
  "A",     2L,  2L,   # extreme underfit floor
  "B",     4L,  4L,   # underfit (old default; oversmooths)
  "C",     8L,  4L,   # approaching adequate
  "D",    13L,  6L,   # the faithful-representation rule
  "E",    16L,  8L,   # current default (rule + margin)
  "F",    24L, 12L,   # M > n_station: overfit onset
  "G",    32L, 16L,   # clear overfit
  "H",    16L,  4L,   # good x, starved y  (same M as I)
  "I",     4L, 16L    # starved x, good y  (same M as H)
) %>% mutate(M = mx * my)

SEEDS <- c(101L, 102L, 103L)

N_CHAINS  <- as.integer(Sys.getenv("SWEEP_CHAINS", "3"))
N_WARMUP  <- as.integer(Sys.getenv("SWEEP_WARMUP", "300"))
N_SAMPLE  <- as.integer(Sys.getenv("SWEEP_SAMPLE", "400"))
PRED_N    <- as.integer(Sys.getenv("SWEEP_PREDN",  "4"))
ADAPT_DELTA <- 0.80
MAX_TREEDEPTH <- 11

# optional subsetting for smoke runs
env_cfg <- Sys.getenv("SWEEP_CONFIGS", ""); env_seed <- Sys.getenv("SWEEP_SEEDS", "")
if (nzchar(env_cfg))  CONFIGS <- CONFIGS %>% filter(label %in% strsplit(env_cfg, ",")[[1]])
if (nzchar(env_seed)) SEEDS   <- as.integer(strsplit(env_seed, ",")[[1]])

OUTPUT_DIR <- "outputs/whale_edna_output_vS1/basis_sweep"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
RESULTS_RDS <- file.path(OUTPUT_DIR, "sweep_results.rds")

cat(sprintf("Sweep: %d configs x %d seeds = %d fits (chains=%d, %d+%d)\n",
            nrow(CONFIGS), length(SEEDS), nrow(CONFIGS) * length(SEEDS),
            N_CHAINS, N_WARMUP, N_SAMPLE))

# -----------------------------------------------------------------------------
# 1. Compile once (basis count is data, so one model serves every config)
# -----------------------------------------------------------------------------
mod <- cmdstan_model("stan/whale_edna_hsgp_vS1.stan")

# resume: load prior results and skip completed (label, seed) pairs
results <- if (file.exists(RESULTS_RDS)) readRDS(RESULTS_RDS) else
  list(diag = tibble(), recovery = tibble(), field = tibble())
done <- if (nrow(results$diag)) results$diag %>% distinct(label, seed) else tibble(label=character(), seed=integer())

# -----------------------------------------------------------------------------
# 2. Metric extraction for one fit
# -----------------------------------------------------------------------------
extract_metrics <- function(fit, sim, lab, mx, my, seed, runtime_s) {
  S  <- sim$meta$n_species
  gp <- sim$truth$gp_params

  # --- sampling health ---
  pd  <- fit$draws(variables = c("mu_sp","gp_sigma","gp_l_raw"))
  ds  <- summarise_draws(pd, rhat = rhat, ess_bulk = ess_bulk, ess_tail = ess_tail)
  dsum <- fit$diagnostic_summary()
  diag <- tibble(label=lab, mx=mx, my=my, M=mx*my, seed=seed,
                 runtime_s = runtime_s,
                 max_rhat  = max(ds$rhat, na.rm=TRUE),
                 min_ess   = min(ds$ess_bulk, na.rm=TRUE),
                 n_divergent = sum(dsum$num_divergent),
                 n_treedepth = sum(dsum$num_max_treedepth),
                 ebfmi_min = min(dsum$ebfmi))

  # --- hyperparameter recovery (lx, ly, gp_sigma, mu) ---
  qsum <- summarise_draws(
    fit$draws(variables = c("gp_l","gp_sigma","mu_sp")),
    mean = mean,
    q025 = ~as.numeric(stats::quantile(.x, 0.025)),
    q975 = ~as.numeric(stats::quantile(.x, 0.975)))
  truth_tbl <- map_dfr(seq_len(S), function(s) tibble(
    species  = s,
    param    = c("lx","ly","gp_sigma","mu"),
    true     = c(gp[[s]]$lx, gp[[s]]$ly, gp[[s]]$sigma, gp[[s]]$mu),
    variable = c(sprintf("gp_l[%d,1]",s), sprintf("gp_l[%d,2]",s),
                 sprintf("gp_sigma[%d]",s), sprintf("mu_sp[%d]",s))))
  rec <- truth_tbl %>%
    left_join(dplyr::select(qsum, variable, mean, q025, q975), by = "variable") %>%
    mutate(label = lab, seed = seed, M = mx*my,
           covered = true >= q025 & true <= q975, bias = mean - true) %>%
    dplyr::select(label, seed, M, species, param, true, mean, q025, q975, covered, bias)

  # --- latent-field recovery: posterior mean log_lambda vs truth ---
  ll   <- summarise_draws(fit$draws(variables = "log_lambda"), mean = mean)
  m    <- stringr::str_match(ll$variable, "\\[(\\d+),(\\d+)\\]")
  post <- matrix(NA_real_, sim$meta$N, S)
  post[cbind(as.integer(m[,2]), as.integer(m[,3]))] <- ll$mean
  truef <- log(sim$truth$lambda_true_si)
  field <- map_dfr(seq_len(S), function(s) {
    tibble(label=lab, seed=seed, M=mx*my, species=s,
           R2   = cor(post[,s], truef[,s])^2,
           rmse = sqrt(mean((post[,s]-truef[,s])^2)),
           bias = mean(post[,s]-truef[,s]))
  })

  list(diag = diag, recovery = rec, field = field)
}

# -----------------------------------------------------------------------------
# 3. Sweep loop (simulate once per seed, refit at each basis count)
# -----------------------------------------------------------------------------
for (seed in SEEDS) {
  sim <- simulate_whale_edna_vS1(seed = seed)
  for (i in seq_len(nrow(CONFIGS))) {
    lab <- CONFIGS$label[i]; mx <- CONFIGS$mx[i]; my <- CONFIGS$my[i]
    if (nrow(done) && any(done$label == lab & done$seed == seed)) {
      cat(sprintf("skip %s / seed %d (done)\n", lab, seed)); next
    }
    cat(sprintf("=== config %s (%dx%d, M=%d) / seed %d ===\n", lab, mx, my, mx*my, seed))
    fd <- format_stan_data_vS1(sim, HSGP_M = c(mx, my), prediction_n = PRED_N)
    sd <- fd$stan_data
    S <- sd$S; D1 <- sd$D1; Mtot <- sd$M
    init_fn <- function() list(
      mu_sp = c(runif(S,1.8,2.0), runif(1,4.0,4.1)),
      gp_sigma = runif(S,0.8,1.2),
      gp_l_raw = matrix(runif(S*D1,0.2,0.3), S, D1, byrow=TRUE),
      z_beta = matrix(0, S, Mtot),
      beta0_phi = 0, gamma0_phi = 2.0, gamma1_phi = 1.0)

    t0 <- Sys.time()
    fit <- tryCatch(
      mod$sample(data = sd, chains = N_CHAINS, parallel_chains = N_CHAINS,
                 iter_warmup = N_WARMUP, iter_sampling = N_SAMPLE,
                 adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH,
                 seed = seed, init = init_fn, refresh = 0, show_messages = FALSE),
      error = function(e) { cat("  FIT ERROR:", conditionMessage(e), "\n"); NULL })
    if (is.null(fit)) next
    runtime_s <- as.numeric(Sys.time() - t0, units = "secs")

    mt <- extract_metrics(fit, sim, lab, mx, my, seed, runtime_s)
    results$diag     <- bind_rows(results$diag,     mt$diag)
    results$recovery <- bind_rows(results$recovery, mt$recovery)
    results$field    <- bind_rows(results$field,    mt$field)
    saveRDS(results, RESULTS_RDS)   # incremental / resumable
    cat(sprintf("  done in %.0fs | maxRhat=%.3f minESS=%.0f div=%d | hake R2=%.3f lx=%.0f\n",
                runtime_s, mt$diag$max_rhat, mt$diag$min_ess, mt$diag$n_divergent,
                mt$field$R2[1], mt$recovery$mean[mt$recovery$species==1 & mt$recovery$param=="lx"]))
  }
}

cat(sprintf("\nSweep complete. %d fits recorded -> %s\n",
            nrow(results$diag), RESULTS_RDS))
