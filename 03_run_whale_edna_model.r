# =============================================================================
# run_whale_edna_model.R
#
# 1. Sources simulate_whale_edna.R to produce synthetic data
# 2. Prepares the Stan data list
#    - GP coordinates use (X, Y, Z_bathy)
#    - Water column effect passed as fixed log_zsample_effect matrix
# 3. Fits whale_edna_hsgp.stan via CmdStanR
# 4. Diagnostics: Rhat, ESS, divergences, energy, traces
# 5. LOO-CV
# 6. Parameter recovery plots
# 7. Posterior predictive checks
# 8. Spatial lambda surfaces (true vs estimated)
# 9. Bathymetric preference recovery
# =============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(patchwork)
library(viridis)

set.seed(42)

# =============================================================================
# 0. Configuration
# =============================================================================

# HSGP basis counts
# X: 300 km domain, min lx ~ 30 km  → 12 basis functions
# Y: 400 km domain, min ly ~ 100 km → 14 basis functions
# Z: 2500 m domain, min lz ~ 150 m  → 12 basis functions
HSGP_M        <- c(6L, 7L, 6L)   # M_total = 2016
HSGP_C        <- c(1.3, 1.3, 1.3)   # boundary factor

N_CHAINS      <- 4
N_WARMUP      <- 250
N_SAMPLE      <- 250
ADAPT_DELTA   <- 0.90
MAX_TREEDEPTH <- 12

OUTPUT_DIR <- "whale_edna_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# =============================================================================
# 1. Simulate data
# =============================================================================
cat("=== Step 1: Simulating data ===\n")
source("01_simulate_whale_edna.R")
sim <- readRDS("whale_edna_sim.rds")

samples       <- sim$design$samples    # X, Y, Z_bathy, Z_sample, depth_idx
stations      <- sim$design$stations   # one row per station
N             <- sim$meta$N
S             <- sim$meta$n_species
R             <- sim$meta$n_qpcr_rep   # 3
K             <- sim$meta$n_mb_rep     # 3
sp_names      <- sim$meta$sp_names
sp_common     <- sim$meta$sp_common
sample_depths <- sim$meta$sample_depths  # 0, 50, 150, 300, 500
vol           <- rep(sim$meta$vol_filtered, N)

# =============================================================================
# 2. Prepare Stan data list
# =============================================================================
cat("=== Step 2: Preparing Stan data ===\n")

# --- GP coordinates: (X, Y, Z_bathy) ---
# Z_bathy spans roughly 10–2500 m across stations
coords_raw <- as.matrix(samples[, c("X", "Y", "Z_bathy")])

# Normalise each dimension to [-1, 1] using fixed physical half-ranges
# X:       0–300 km  → centre 150,  half-range 150
# Y:       0–400 km  → centre 200,  half-range 200
# Z_bathy: 0–2500 m  → centre 1250, half-range 1250
coord_center <- c(150,  200,  1250)
coord_scale  <- c(150,  200,  1250)

coords_norm <- sweep(coords_raw, 2, coord_center, "-")
coords_norm <- sweep(coords_norm, 2, coord_scale,  "/")

L_hsgp  <- HSGP_C
M_total <- prod(HSGP_M)   # 2016

cat(sprintf("  Coordinate ranges after normalisation:\n"))
cat(sprintf("    X:       [%.2f, %.2f]\n", min(coords_norm[,1]), max(coords_norm[,1])))
cat(sprintf("    Y:       [%.2f, %.2f]\n", min(coords_norm[,2]), max(coords_norm[,2])))
cat(sprintf("    Z_bathy: [%.2f, %.2f]\n", min(coords_norm[,3]), max(coords_norm[,3])))

# --- Water column effect: log_zsample_effect ---
# Pre-computed from simulation truth. In a real analysis this would either
# be fixed from literature values or estimated as parameters in Stan.
# Passed as an N × S matrix of log-scale multipliers.
log_zsample_effect <- log(sim$truth$zsample_effect)   # N × S

cat(sprintf("  log_zsample_effect range: [%.2f, %.2f]\n",
            min(log_zsample_effect), max(log_zsample_effect)))

# --- qPCR data (hake only, s=1) ---
n_detect    <- sim$observed$qpcr_n_detect     # N integer
any_detect  <- sim$observed$qpcr_any_detect   # N {0,1}
mean_ct_obs <- sim$observed$qpcr_mean_ct      # N real
mean_ct_obs[is.na(mean_ct_obs)] <- 0.0        # Stan ignores when any_detect==0

# --- Metabarcoding (all species, K replicates) ---
mb_reads_arr <- sim$observed$mb_reads_rep    # N × S × K array
mb_total_arr <- sim$observed$mb_total_rep    # N × K matrix

# Guard volume — replace any zeros with a small positive value
vol <- rep(sim$meta$vol_filtered, N)
vol <- pmax(vol, 1e-6)
# Guard log_zsample_effect — replace any -Inf or extreme values
log_zsample_effect <- log(sim$truth$zsample_effect)
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect < 0] <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect > 0] <-  10.0
log_zsample_effect[is.nan(log_zsample_effect)] <- -10.0
cat("vol range after guard:", range(vol), "\n")
cat("log_zsample_effect range after guard:", range(log_zsample_effect), "\n")

stan_data <- list(
  N                  = N,
  S                  = S,
  R                  = R,
  K                  = K,
  M                  = M_total,
  coords             = coords_norm,
  coord_scale        = coord_scale,
  L_hsgp             = L_hsgp,
  m_hsgp             = HSGP_M,
  vol                = vol,
  log_zsample_effect = log_zsample_effect,
  # qPCR
  n_detect           = n_detect,
  mean_ct_obs        = mean_ct_obs,
  any_detect         = any_detect,
  # Metabarcoding
  mb_reads           = mb_reads_arr,
  mb_total           = mb_total_arr
)

cat(sprintf("  Stan data prepared: N=%d, S=%d, R=%d, K=%d, M=%d\n",
            N, S, R, K, M_total))

# =============================================================================
# 3. Compile and fit
# =============================================================================
cat("=== Step 3: Compiling Stan model ===\n")
mod <- cmdstan_model("whale_edna_hsgp.stan")

cat("=== Step 3: Fitting model ===\n")
init_fn <- function() {
  list(
    gp_sigma = rep(0.5, S),
    gp_l     = matrix(c(50, 150, 300), nrow = S, ncol = 3, byrow = TRUE),
    z_beta   = matrix(0, nrow = S, ncol = M_total),
    mu_sp    = rep(2.0, S),          # log(~7 copies/L) — reasonable start
    kappa    = 0.8,
    alpha_ct = 38.0,
    beta_ct  = 3.32,
    sigma_ct = 0.5,
    p_zi     = 0.3,
    phi_bb   = rep(10.0, S)
  )
}
fit <- mod$sample(
  data            = stan_data,
  chains          = N_CHAINS,
  parallel_chains = N_CHAINS,
  iter_warmup     = N_WARMUP,
  iter_sampling   = N_SAMPLE,
  adapt_delta     = ADAPT_DELTA,
  max_treedepth   = MAX_TREEDEPTH,
  seed            = 42,
  init            = init_fn,          # <-- add this line
  output_dir      = OUTPUT_DIR,
  show_messages   = TRUE
)

# =============================================================================
# 4. Diagnostics
# =============================================================================
cat("=== Step 4: Diagnostics ===\n")

fit$cmdstan_diagnose()

draws_arr <- fit$draws(format = "array")
draws_df  <- fit$draws(format = "df")

diag_summary <- summarise_draws(
  fit$draws(),
  mean, sd,
  ~quantile(.x, probs = c(0.025, 0.975)),
  Rhat     = rhat,
  ESS_bulk = ess_bulk,
  ESS_tail = ess_tail
)

rhat_bad <- diag_summary %>% filter(Rhat > 1.01)
ess_bad  <- diag_summary %>% filter(ESS_bulk < 400 | ESS_tail < 400)

cat(sprintf("  Parameters with Rhat > 1.01 : %d\n", nrow(rhat_bad)))
cat(sprintf("  Parameters with ESS < 400   : %d\n", nrow(ess_bad)))
if (nrow(rhat_bad) > 0) print(rhat_bad %>% select(variable, Rhat))

sampler_diag <- fit$sampler_diagnostics(format = "df")
n_div <- sum(sampler_diag$divergent__)
n_max <- sum(sampler_diag$treedepth__ >= MAX_TREEDEPTH)
cat(sprintf("  Divergences     : %d\n", n_div))
cat(sprintf("  Max treedepth   : %d\n", n_max))

# Energy plot
energy_plot <- mcmc_nuts_energy(sampler_diag) +
  ggtitle("NUTS Energy — Bayesian Fraction of Missing Information")
ggsave(file.path(OUTPUT_DIR, "diag_energy.png"),
       energy_plot, width = 8, height = 4)

# Trace plots for key scalar parameters
scalar_pars <- c(
  paste0("gp_sigma[", 1:S, "]"),
  paste0("mu_sp[",    1:S, "]"),
  paste0("kappa"),
  paste0("p_zi"),
  paste0("phi_bb[",   1:S, "]")
)
trace_p <- mcmc_trace(draws_arr, pars = scalar_pars) +
  ggtitle("Trace plots — scalar parameters")
ggsave(file.path(OUTPUT_DIR, "diag_trace.png"),
       trace_p, width = 12, height = 10)

# =============================================================================
# 5. LOO-CV
# =============================================================================
cat("=== Step 5: LOO-CV ===\n")

log_lik_qpcr <- fit$draws("log_lik_qpcr", format = "matrix")
log_lik_mb   <- fit$draws("log_lik_mb",   format = "matrix")

loo_qpcr <- loo(log_lik_qpcr, cores = 2)
loo_mb   <- loo(log_lik_mb,   cores = 2)

cat("  LOO — qPCR:\n");         print(loo_qpcr)
cat("  LOO — Metabarcoding:\n"); print(loo_mb)

png(file.path(OUTPUT_DIR, "loo_pareto_k.png"), width = 900, height = 400)
par(mfrow = c(1, 2))
plot(loo_qpcr, main = "LOO Pareto-k: qPCR (hake)")
plot(loo_mb,   main = "LOO Pareto-k: Metabarcoding")
dev.off()

# =============================================================================
# 6. Parameter recovery
# =============================================================================
cat("=== Step 6: Parameter recovery ===\n")

scalar_summary <- diag_summary %>%
  filter(variable %in% scalar_pars) %>%
  mutate(
    param  = str_extract(variable, "^[a-z_]+"),
    sp_idx = as.integer(str_extract(variable, "[0-9]+"))
  )

# True values
true_vals <- tibble(
  param  = rep(c("gp_sigma", "mu_sp", "phi_bb"), each = S),
  sp_idx = rep(1:S, 3),
  true   = c(
    sapply(sim$truth$gp_params, `[[`, "sigma"),
    sapply(sim$truth$gp_params, `[[`, "mu"),
    sim$truth$phi_bb
  )
) %>%
  bind_rows(
    tibble(param = "kappa", sp_idx = 1L,
           true  = sim$truth$qpcr_params$kappa),
    tibble(param = "p_zi",  sp_idx = 1L,
           true  = sim$truth$p_zi[1])
  )

recovery_df <- scalar_summary %>%
  left_join(true_vals, by = c("param", "sp_idx")) %>%
  filter(!is.na(true)) %>%
  mutate(sp_name = ifelse(is.na(sp_idx), "hake",
                          sp_common[sp_idx]))

p_recovery <- ggplot(recovery_df,
                     aes(x = true, y = mean,
                         ymin = `2.5%`, ymax = `97.5%`,
                         colour = sp_name, shape = sp_name)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50") +
  geom_errorbar(width = 0, alpha = 0.6) +
  geom_point(size = 3) +
  facet_wrap(~param, scales = "free") +
  scale_colour_viridis_d(option = "D", end = 0.85) +
  labs(
    title  = "Parameter recovery: true vs posterior mean (95% CI)",
    x = "True value", y = "Posterior mean ± 95% CI",
    colour = "Species", shape = "Species"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUTPUT_DIR, "param_recovery.png"),
       p_recovery, width = 12, height = 8)

# Length-scale recovery
ls_pars <- paste0("gp_l[", rep(1:S, each=3), ",", rep(1:3, S), "]")

ls_summary <- diag_summary %>%
  filter(variable %in% ls_pars) %>%
  mutate(
    sp_idx   = as.integer(str_extract(variable, "(?<=\\[)[0-9]+")),
    dim_idx  = as.integer(str_extract(variable, "(?<=,)[0-9](?=\\])")),
    dim_name = c("lx (km)", "ly (km)", "lz_bathy (m)")[dim_idx],
    sp_name  = sp_common[sp_idx]
  )

true_ls <- bind_rows(lapply(seq_len(S), function(s) {
  p <- sim$truth$gp_params[[s]]
  tibble(
    sp_idx   = s,
    dim_idx  = 1:3,
    true     = c(p$lx, p$ly, p$lz),
    dim_name = c("lx (km)", "ly (km)", "lz_bathy (m)"),
    sp_name  = sp_common[s]
  )
}))

ls_df <- ls_summary %>%
  left_join(true_ls, by = c("sp_idx", "dim_idx", "dim_name", "sp_name"))

p_ls <- ggplot(ls_df,
               aes(x = true, y = mean,
                   ymin = `2.5%`, ymax = `97.5%`,
                   colour = sp_name, shape = dim_name)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50") +
  geom_errorbar(width = 0, alpha = 0.6) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "C", end = 0.85) +
  labs(
    title  = "Length-scale recovery (lx, ly in km; lz_bathy in m)",
    x = "True length-scale", y = "Posterior mean ± 95% CI",
    colour = "Species", shape = "Dimension"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUTPUT_DIR, "lengthscale_recovery.png"),
       p_ls, width = 8, height = 5)

# =============================================================================
# 7. Posterior predictive checks
# =============================================================================
cat("=== Step 7: Posterior predictive checks ===\n")

pp_n_detect <- fit$draws("pp_n_detect", format = "matrix")
pp_mean_ct  <- fit$draws("pp_mean_ct",  format = "matrix")
pp_mb_reads <- fit$draws("pp_mb_reads", format = "matrix")

# Reshape draws × (N*S) → list of draws × N per species
# Stan stores [N, S] column-major: index = (s-1)*N + i
reshape_ns <- function(mat, N, S) {
  lapply(1:S, function(s) mat[, ((s-1)*N + 1):(s*N)])
}

# mb_reads is [N, S, K]: Stan column-major over all three dims
# Reshape to list[s][k] of draws × N
reshape_nsk <- function(mat, N, S, K) {
  lapply(1:S, function(s) {
    lapply(1:K, function(k) {
      col_start <- ((k-1)*S*N) + ((s-1)*N) + 1
      mat[, col_start:(col_start + N - 1)]
    })
  })
}

pp_det_list <- list(pp_n_detect)   # N-length, hake only
pp_ct_list  <- list(pp_mean_ct)
pp_mb_list  <- reshape_nsk(pp_mb_reads, N, S, K)

ppc_plots <- list()

# qPCR (hake only)
y_det <- n_detect
ppc_plots[["det_hake"]] <-
  ppc_bars(y_det, pp_det_list[[1]][1:200, ]) +
  ggtitle("PPC: qPCR n_detect — Pacific hake") +
  theme_bw(base_size = 11)

idx_det <- which(any_detect == 1)
if (length(idx_det) > 5) {
  y_ct <- mean_ct_obs[idx_det]
  ppc_plots[["ct_hake"]] <-
    ppc_dens_overlay(y_ct, pp_ct_list[[1]][1:200, idx_det]) +
    ggtitle("PPC: mean Ct (detected only) — Pacific hake") +
    theme_bw(base_size = 11)
}

# Metabarcoding (all species, replicate 1 for display)
for (s in seq_len(S)) {
  y_mb <- mb_reads_arr[, s, 1]
  ppc_plots[[paste0("mb_", s)]] <-
    ppc_dens_overlay(y_mb, pp_mb_list[[s]][[1]][1:100, ]) +
    ggtitle(sprintf("PPC: MB reads rep 1 — %s", sp_common[s])) +
    theme_bw(base_size = 11)
}

ppc_grid <- wrap_plots(ppc_plots, ncol = 2)
ggsave(file.path(OUTPUT_DIR, "ppc_all.png"),
       ppc_grid, width = 12,
       height = 4 * ceiling(length(ppc_plots) / 2))

# =============================================================================
# 8. Spatial lambda surface: true vs estimated
# =============================================================================
cat("=== Step 8: Spatial lambda plots ===\n")

lambda_hat_draws <- fit$draws("lambda_hat", format = "matrix")
lambda_hat_list  <- reshape_ns(lambda_hat_draws, N, S)

lambda_post_mean <- matrix(NA, N, S)
lambda_post_sd   <- matrix(NA, N, S)
for (s in seq_len(S)) {
  lambda_post_mean[, s] <- colMeans(lambda_hat_list[[s]])
  lambda_post_sd[, s]   <- apply(lambda_hat_list[[s]], 2, sd)
}

# Plot at surface samples only (Z_sample == 0)
surf_idx <- which(samples$Z_sample == 0)

plot_df <- bind_rows(lapply(seq_len(S), function(s) {
  data.frame(
    X           = samples$X[surf_idx],
    Y           = samples$Y[surf_idx],
    Z_bathy     = samples$Z_bathy[surf_idx],
    species     = sp_common[s],
    lambda_true = lambda_true[surf_idx, s],
    lambda_est  = lambda_post_mean[surf_idx, s],
    lambda_sd   = lambda_post_sd[surf_idx, s]
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

p_scatter <- ggplot(plot_df,
                    aes(x = log(lambda_true + 0.01),
                        y = log(lambda_est  + 0.01),
                        colour = species)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.4, size = 1.2) +
  facet_wrap(~species) +
  scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
  labs(
    title = "True vs posterior mean log-lambda (surface samples)",
    x = "log(true λ)", y = "log(posterior mean λ)"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUTPUT_DIR, "lambda_scatter.png"),
       p_scatter, width = 12, height = 4)
# =============================================================================
# Step 9 (continued): Bathymetric preference recovery
# =============================================================================

p_bathy_rec <- ggplot(bathy_recovery_df,
                      aes(x = Z_bathy)) +
  geom_point(aes(y = log(lambda_true + 0.01)),
             colour = "grey60", size = 1.0, alpha = 0.4) +
  geom_point(aes(y = log(lambda_est + 0.01), colour = species),
             size = 1.2, alpha = 0.6) +
  facet_wrap(~species) +
  scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
  labs(
    title    = "Bathymetric preference: true (grey) vs estimated (colour)",
    subtitle = "log(λ) vs Z_bathy at surface samples; grey = true, colour = posterior mean",
    x = "Bottom depth Z_bathy (m)",
    y = "log(λ)"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUTPUT_DIR, "bathy_preference_recovery.png"),
       p_bathy_rec, width = 12, height = 4)

# Smooth comparison: true vs estimated mean function over Z_bathy gradient
zbathy_seq <- seq(10, 2500, by = 20)

bathy_fn_df <- bind_rows(lapply(seq_len(S), function(s) {
  p <- sim$truth$gp_params[[s]]
  
  zbathy_pref_fn <- function(Z_bathy, mu_z, sd_z, amp) {
    raw <- dnorm(Z_bathy, mu_z, sd_z)
    raw <- raw / max(raw)
    amp * (raw - 0.5)
  }
  
  z_mu_true <- zbathy_pref_fn(zbathy_seq,
                              p$zbathy_pref_mu,
                              p$zbathy_pref_sd,
                              p$zbathy_pref_amp)
  data.frame(
    Z_bathy  = zbathy_seq,
    log_lambda = p$mu + z_mu_true,
    type     = "True mean function",
    species  = sp_common[s]
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

p_bathy_fn <- ggplot() +
  geom_line(
    data      = bathy_fn_df,
    aes(x = Z_bathy, y = log_lambda, colour = species),
    linewidth = 1.2,
    linetype  = "solid"
  ) +
  geom_vline(
    xintercept = 200,
    linetype   = "dashed",
    colour     = "grey50",
    linewidth  = 0.5
  ) +
  annotate(
    "text", x = 210, y = max(bathy_fn_df$log_lambda) * 0.98,
    label    = "shelf break\n(~200 m)",
    colour   = "grey40",
    size     = 3.0,
    hjust    = 0,
    fontface = "italic"
  ) +
  scale_colour_manual(values = sp_colours, name = NULL,
                      labels = sp_common) +
  scale_x_continuous(
    breaks = c(0, 200, 500, 1000, 1500, 2000, 2500),
    expand = c(0.02, 0)
  ) +
  labs(
    title    = "True bathymetric mean function by species",
    subtitle = "log(λ) = mu_s + zbathy_pref(Z_bathy)",
    x = "Bottom depth Z_bathy (m)",
    y = "log(λ)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text     = element_text(face = "italic", size = 9)
  )

ggsave(file.path(OUTPUT_DIR, "bathy_mean_function.png"),
       p_bathy_fn, width = 9, height = 5)

# =============================================================================
# 10. Save diagnostics summary to CSV
# =============================================================================
cat("=== Step 10: Saving diagnostics ===\n")

write_csv(diag_summary,
          file.path(OUTPUT_DIR, "diagnostics_summary.csv"))

if (nrow(rhat_bad) > 0) {
  write_csv(rhat_bad,
            file.path(OUTPUT_DIR, "rhat_flagged.csv"))
  cat(sprintf("  WARNING: %d parameters with Rhat > 1.01 written to rhat_flagged.csv\n",
              nrow(rhat_bad)))
}

if (nrow(ess_bad) > 0) {
  write_csv(ess_bad,
            file.path(OUTPUT_DIR, "ess_flagged.csv"))
  cat(sprintf("  WARNING: %d parameters with ESS < 400 written to ess_flagged.csv\n",
              nrow(ess_bad)))
}

# Sampler diagnostics summary
sampler_summary <- data.frame(
  metric = c("Divergences", "Max treedepth hits",
             "Mean accept stat", "Mean stepsize"),
  value  = c(
    n_div,
    n_max,
    mean(sampler_diag$accept_stat__),
    mean(sampler_diag$stepsize__)
  )
)
write_csv(sampler_summary,
          file.path(OUTPUT_DIR, "sampler_summary.csv"))
print(sampler_summary)

# =============================================================================
# 11. Session info
# =============================================================================
cat("=== Step 11: Session info ===\n")

sink(file.path(OUTPUT_DIR, "session_info.txt"))
print(sessionInfo())
sink()

cat("\n=== All done. Outputs written to:", OUTPUT_DIR, "===\n")
cat("Files produced:\n")
for (f in list.files(OUTPUT_DIR)) {
  cat(sprintf("  %s\n", f))
}
