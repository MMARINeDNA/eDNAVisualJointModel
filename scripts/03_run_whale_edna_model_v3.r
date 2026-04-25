# =============================================================================
# run_whale_edna_model.R
#
# Runner script for whale_edna_hsgp.stan (v2).
# Reads whale_edna_sim.rds produced by simulate_whale_edna.R (v1).
#
# Spatial structure is handled exclusively by the 3-D anisotropic HSGP
# over (X, Y, Z_bathy). No additive zbathy or spatial mean functions.
#
# Data preparation:
#   - GP coordinates: (X, Y, Z_bathy), normalised to [-1, 1]
#   - log_zsample_effect: fixed N × S water-column offset, pre-computed in R
#   - qPCR: long-form N*R rows (one per replicate), hake only
#   - MB: long-form N*K rows (one per replicate), filtered to rows where
#         at least one species has reads > 0
# =============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(patchwork)
library(viridis)

set.seed(42)

# Need to feed the model fixed values of alpha and beta in this version
# in version with standard curve estimation those values will be estimated

# =============================================================================
# 0. Configuration
# =============================================================================

# HSGP basis terms per dimension.
# lx_min ~ 30 km in 300 km domain  → ~10 terms
# ly_min ~ 100 km in 400 km domain → ~8 terms
# lz_min ~ 100 m in 2500 m domain  → ~8 terms
HSGP_M        <- c(10L, 8L, 8L)     # M_total = 640
HSGP_C        <- c(1.5, 1.5, 1.5)   # boundary extension (>= 1.5 recommended)

N_CHAINS      <- 4
N_WARMUP      <- 500
N_SAMPLE      <- 500
ADAPT_DELTA   <- 0.90
MAX_TREEDEPTH <- 12

OUTPUT_DIR <- "whale_edna_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# =============================================================================
# 1. Load simulation
# =============================================================================
cat("=== Step 1: Loading simulation ===\n")
sim <- readRDS("whale_edna_sim.rds")

samples         <- sim$design$samples
stations        <- sim$design$stations
gp_params       <- sim$truth$gp_params
lambda_true_si  <- sim$truth$lambda_true_si    # N × S
zsample_effect  <- sim$truth$zsample_effect    # N × S
C_obs_si        <- sim$truth$C_obs_si          # N × S
qpcr_detect_raw <- sim$observed$qpcr_detect_raw  # N*R × S (col 1 = hake)
qpcr_ct_raw     <- sim$observed$qpcr_ct_raw      # N*R × S
mb_reads_rep    <- sim$observed$mb_reads_rep     # S × N*R (from rmultinom)

N            <- sim$meta$N
S            <- sim$meta$n_species
R            <- sim$meta$n_qpcr_rep    # 3
K            <- sim$meta$n_mb_rep      # 3
sp_common    <- sim$meta$sp_common
vol_aliquot  <- 2                       # microlitres (from simulation script)
conv_factor  <- 10                      # animals → copies

cat(sprintf("  Loaded: N=%d samples, S=%d species, R=%d qPCR reps, K=%d MB reps\n",
            N, S, R, K))

# =============================================================================
# 2. Water column log-effect (fixed offset, pre-computed from simulation)
# =============================================================================
cat("=== Step 2: Preparing log_zsample_effect ===\n")

log_zsample_effect <- log(as.matrix(zsample_effect))   # N × S

# Guard against -Inf / NaN at extreme depths
log_zsample_effect[is.nan(log_zsample_effect)]                           <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect < 0] <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect > 0] <-  10.0

cat(sprintf("  log_zsample_effect range: [%.2f, %.2f]\n",
            min(log_zsample_effect), max(log_zsample_effect)))
stopifnot(all(is.finite(log_zsample_effect)))

# =============================================================================
# 3. GP coordinates — normalise to [-1, 1]
# =============================================================================
cat("=== Step 3: Normalising coordinates ===\n")

coords_raw <- cbind(
  as.numeric(samples[["X"]]),
  as.numeric(samples[["Y"]]),
  as.numeric(samples[["Z_bathy"]])
)

# Fixed physical ranges (interpretable; consistent across subsets)
# X:       0 – 300 km  → centre = 150, half-range = 150
# Y:       0 – 400 km  → centre = 200, half-range = 200
# Z_bathy: 0 – 2500 m  → centre = 1250, half-range = 1250
coord_centre <- c(150,  200,  1250)
coord_scale  <- c(150,  200,  1250)

coords_norm <- sweep(coords_raw, 2, coord_centre, "-")
coords_norm <- sweep(coords_norm, 2, coord_scale,  "/")

cat(sprintf("  X_norm:       [%.3f, %.3f]\n", min(coords_norm[,1]), max(coords_norm[,1])))
cat(sprintf("  Y_norm:       [%.3f, %.3f]\n", min(coords_norm[,2]), max(coords_norm[,2])))
cat(sprintf("  Z_bathy_norm: [%.3f, %.3f]\n", min(coords_norm[,3]), max(coords_norm[,3])))
stopifnot(all(is.finite(coords_norm)))

# =============================================================================
# 4. qPCR data — long form (N * R rows)
# =============================================================================
cat("=== Step 4: Preparing qPCR data ===\n")

N_qpcr_long     <- N * R
qpcr_sample_idx <- rep(seq_len(N), each = R)

qpcr_detect_vec <- as.integer(qpcr_detect_raw[, 1])
qpcr_ct_vec     <- as.numeric(qpcr_ct_raw[, 1])
qpcr_ct_vec[is.na(qpcr_ct_vec)] <- 0.0   # Stan requires numeric; ignored when detect == 0

stopifnot(all(qpcr_detect_vec %in% c(0L, 1L)))
stopifnot(length(qpcr_detect_vec) == N_qpcr_long)

cat(sprintf("  N_qpcr_long: %d\n", N_qpcr_long))
cat(sprintf("  qPCR detection rate (hake): %.1f%%\n",
            100 * mean(qpcr_detect_vec)))

# =============================================================================
# 5. Metabarcoding data — long form (N * K rows), zero-read rows filtered
# =============================================================================
cat("=== Step 5: Preparing metabarcoding data ===\n")

# mb_reads_rep from rmultinom is S × (N*K); transpose to (N*K) × S
mb_reads_long_full <- t(mb_reads_rep)
storage.mode(mb_reads_long_full) <- "integer"

mb_sample_idx_full <- rep(seq_len(N), each = K)
mb_total_full      <- rowSums(mb_reads_long_full)

# Filter rows where ALL species have zero reads
keep_mb      <- mb_total_full > 0
n_dropped_mb <- sum(!keep_mb)

mb_reads_long <- mb_reads_long_full[keep_mb, , drop = FALSE]
mb_sample_idx <- mb_sample_idx_full[keep_mb]
mb_total_vec  <- mb_total_full[keep_mb]
N_mb_long     <- nrow(mb_reads_long)
storage.mode(mb_total_vec) <- "integer"

cat(sprintf("  MB rows (N*K = %d): %d kept, %d dropped (zero reads)\n",
            N * K, N_mb_long, n_dropped_mb))
cat(sprintf("  MB read depth range: %d – %d\n",
            min(mb_total_vec), max(mb_total_vec)))

stopifnot(all(mb_total_vec > 0))
stopifnot(all(mb_reads_long >= 0L))
stopifnot(ncol(mb_reads_long) == S)

# =============================================================================
# 6. Assemble Stan data list
# =============================================================================
cat("=== Step 6: Assembling Stan data list ===\n")

M_total <- prod(HSGP_M)

stan_data <- list(
  
  # Dimensions
  N  = N,
  S  = S,
  M  = M_total,
  
  # Coordinates (normalised to [-1, 1])
  coords      = coords_norm,
  coord_scale = coord_scale,
  
  # HSGP
  L_hsgp = HSGP_C,
  m_hsgp = HSGP_M,
  
  # Water column offset (fixed)
  log_zsample_effect = log_zsample_effect,
  
  # Conversion factor
  log_conv_factor = log(conv_factor),
  
  # Aliquot volume (microlitres)
  vol_aliquot = vol_aliquot,
  
  # qPCR (hake only, long form)
  N_qpcr_long     = N_qpcr_long,
  qpcr_sample_idx = qpcr_sample_idx,
  qpcr_detect     = qpcr_detect_vec,
  qpcr_ct         = qpcr_ct_vec,
  
  # Metabarcoding (all species, long form, zero-read rows removed)
  N_mb_long     = N_mb_long,
  mb_sample_idx = mb_sample_idx,
  mb_reads      = mb_reads_long,
  mb_total      = mb_total_vec,
  
  # Zero-inflation switch
  use_zi = 1L,
  
  # Prior hyperparameters
  prior_mu_sp_mu     =  2.0,
  prior_mu_sp_sig    =  1.5,
  prior_gp_sigma_sig =  1.5,
  prior_gp_lx_mu     = 50.0,  prior_gp_lx_sig  = 40.0,
  prior_gp_ly_mu     = 150.0, prior_gp_ly_sig  = 80.0,
  prior_gp_lz_mu     = 300.0, prior_gp_lz_sig  = 150.0
)

# Sanity checks
stopifnot(nrow(coords_norm)        == N)
stopifnot(nrow(log_zsample_effect) == N)
stopifnot(ncol(log_zsample_effect) == S)
stopifnot(length(qpcr_sample_idx)  == N_qpcr_long)
stopifnot(length(mb_sample_idx)    == N_mb_long)
stopifnot(nrow(mb_reads_long)      == N_mb_long)

cat(sprintf("  Stan data ready: N=%d  S=%d  M=%d  N_qpcr=%d  N_mb=%d\n",
            N, S, M_total, N_qpcr_long, N_mb_long))

# =============================================================================
# 7. Compile and fit
# =============================================================================
cat("=== Step 7: Compiling Stan model ===\n")

mod <- cmdstan_model(
  "whale_edna_hsgp_V2.stan",
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
    alpha_ct   = 38.0,
    beta_ct    = 1.44,
    sigma_ct   = 0.50,
    beta0_phi  = rep(2.0, S),
    gamma0_phi = rep(5.0, S),
    gamma1_phi = rep(1.0, S)
  )
}

cat("=== Step 7: Fitting model ===\n")
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
cat("Fit saved.\n")
# =============================================================================
# 8. Diagnostics
# =============================================================================
cat("=== Step 8: Diagnostics ===\n")

fit$cmdstan_diagnose()

draws_arr <- fit$draws(format = "array")
diag_sum  <- summarise_draws(
  fit$draws(),
  mean, sd,
  ~quantile(.x, probs = c(0.025, 0.975)),
  Rhat     = rhat,
  ESS_bulk = ess_bulk,
  ESS_tail = ess_tail
)

rhat_bad <- diag_sum %>% filter(Rhat > 1.01)
ess_bad  <- diag_sum %>% filter(ESS_bulk < 400 | ESS_tail < 400)
cat(sprintf("  Rhat > 1.01 : %d parameters\n", nrow(rhat_bad)))
cat(sprintf("  ESS < 400   : %d parameters\n", nrow(ess_bad)))
if (nrow(rhat_bad) > 0) print(rhat_bad %>% select(variable, Rhat))

sampler_diag <- fit$sampler_diagnostics(format = "df")
cat(sprintf("  Divergences       : %d\n", sum(sampler_diag$divergent__)))
cat(sprintf("  Max treedepth hits: %d\n",
            sum(sampler_diag$treedepth__ >= MAX_TREEDEPTH)))

# Energy plot
energy_plot <- mcmc_nuts_energy(sampler_diag) +
  ggtitle("NUTS Energy")
ggsave(file.path(OUTPUT_DIR, "diag_energy.png"),
       energy_plot, width = 8, height = 4)

# Scalar parameter trace plots
scalar_pars <- c(
  paste0("mu_sp[",    1:S, "]"),
  paste0("gp_sigma[", 1:S, "]"),
  "kappa", "alpha_ct", "beta_ct", "sigma_ct",
  paste0("beta0_phi[",  1:S, "]"),
  paste0("gamma0_phi[", 1:S, "]"),
  paste0("gamma1_phi[", 1:S, "]")
)

trace_plot <- mcmc_trace(draws_arr, pars = scalar_pars) +
  ggtitle("Trace plots")
ggsave(file.path(OUTPUT_DIR, "diag_trace.png"),
       trace_plot, width = 14, height = 10)

print(fit$summary(scalar_pars))

# =============================================================================
# 9. LOO-CV
# =============================================================================
cat("=== Step 9: LOO-CV ===\n")

loo_qpcr <- loo(fit$draws("log_lik_qpcr", format = "matrix"), cores = 2)
loo_mb   <- loo(fit$draws("log_lik_mb",   format = "matrix"), cores = 2)
cat("qPCR LOO:\n");          print(loo_qpcr)
cat("Metabarcoding LOO:\n"); print(loo_mb)

# =============================================================================
# 10. Parameter recovery
# =============================================================================
cat("=== Step 10: Parameter recovery ===\n")

true_scalar <- tibble(
  variable = c(
    paste0("mu_sp[",    1:S, "]"),
    paste0("gp_sigma[", 1:S, "]"),
    "kappa", "alpha_ct", "beta_ct", "sigma_ct"
  ),
  true = c(
    sapply(gp_params, `[[`, "mu"),
    sapply(gp_params, `[[`, "sigma"),
    sim$truth$qpcr_params$kappa,
    sim$truth$qpcr_params$alpha_ct,
    sim$truth$qpcr_params$beta_ct,
    sim$truth$qpcr_params$sigma_ct
  )
)

recovery_df <- diag_sum %>%
  filter(variable %in% true_scalar$variable) %>%
  left_join(true_scalar, by = "variable")

p_recovery <- ggplot(recovery_df,
                     aes(x = true, y = mean,
                         ymin = `2.5%`, ymax = `97.5%`)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50") +
  geom_errorbar(width = 0, alpha = 0.6) +
  geom_point(size = 3, colour = "#3B528B") +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Parameter recovery: true vs posterior mean ± 95% CI",
       x = "True value", y = "Posterior mean") +
  theme_bw(base_size = 11)

ggsave(file.path(OUTPUT_DIR, "param_recovery.png"),
       p_recovery, width = 14, height = 8)

# Length-scale recovery
ls_pars <- paste0("gp_l[", rep(1:S, each = 3), ",", rep(1:3, S), "]")

true_ls <- bind_rows(lapply(seq_len(S), function(s) {
  p <- gp_params[[s]]
  tibble(
    variable = paste0("gp_l[", s, ",", 1:3, "]"),
    true     = c(p$lx, p$ly, p$lz),
    dim_name = c("lx (km)", "ly (km)", "lz_bathy (m)"),
    sp_name  = sp_common[s]
  )
}))

ls_recovery_df <- diag_sum %>%
  filter(variable %in% ls_pars) %>%
  left_join(true_ls, by = "variable")

p_ls <- ggplot(ls_recovery_df,
               aes(x = true, y = mean,
                   ymin = `2.5%`, ymax = `97.5%`,
                   colour = sp_name, shape = dim_name)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50") +
  geom_errorbar(width = 0, alpha = 0.6) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "D", end = 0.85) +
  labs(title  = "Length-scale recovery",
       x      = "True length-scale",
       y      = "Posterior mean ± 95% CI",
       colour = "Species", shape = "Dimension") +
  theme_bw(base_size = 11)

ggsave(file.path(OUTPUT_DIR, "lengthscale_recovery.png"),
       p_ls, width = 9, height = 5)

# =============================================================================
# 11. Posterior predictive checks
# =============================================================================
cat("=== Step 11: Posterior predictive checks ===\n")

pp_qpcr_detect <- fit$draws("pp_qpcr_detect", format = "matrix")
pp_qpcr_ct     <- fit$draws("pp_qpcr_ct",     format = "matrix")
pp_mb_reads    <- fit$draws("pp_mb_reads",    format = "matrix")

# Stan stores array[N_mb_long, S] column-major: species s occupies
# columns ((s-1)*N_mb_long + 1) : (s*N_mb_long)
extract_mb_species <- function(mat, s, N_mb) {
  mat[, ((s - 1) * N_mb + 1):(s * N_mb), drop = FALSE]
}

ppc_plots <- list()

# qPCR detection (0/1)
ppc_plots[["qpcr_detect"]] <-
  ppc_bars(
    y    = qpcr_detect_vec,
    yrep = pp_qpcr_detect[1:min(200, nrow(pp_qpcr_detect)), ]
  ) +
  ggtitle("PPC: qPCR detection (0/1) — Pacific hake") +
  theme_bw(base_size = 11)

# Ct values (detected replicates only)
det_idx <- which(qpcr_detect_vec == 1)
if (length(det_idx) >= 5) {
  ppc_plots[["qpcr_ct"]] <-
    ppc_dens_overlay(
      y    = qpcr_ct_vec[det_idx],
      yrep = pp_qpcr_ct[1:min(200, nrow(pp_qpcr_ct)), det_idx, drop = FALSE]
    ) +
    ggtitle("PPC: Ct values (detected only) — Pacific hake") +
    theme_bw(base_size = 11)
}

# MB read counts per species
for (s in seq_len(S)) {
  pp_mb_s <- extract_mb_species(pp_mb_reads, s, N_mb_long)
  ppc_plots[[paste0("mb_sp", s)]] <-
    ppc_dens_overlay(
      y    = mb_reads_long[, s],
      yrep = pp_mb_s[1:min(100, nrow(pp_mb_s)), , drop = FALSE]
    ) +
    ggtitle(sprintf("PPC: MB read counts — %s", sp_common[s])) +
    theme_bw(base_size = 11)
}

ggsave(file.path(OUTPUT_DIR, "ppc_all.png"),
       wrap_plots(ppc_plots, ncol = 2),
       width  = 12,
       height = 4 * ceiling(length(ppc_plots) / 2))

# =============================================================================
# 12. Spatial lambda: true vs posterior mean
# =============================================================================
cat("=== Step 12: Spatial lambda plots ===\n")

# Extract posterior draws of lambda_hat — Stan stores N × S column-major
# so species s occupies columns ((s-1)*N + 1) : (s*N)
lambda_hat_draws <- fit$draws("lambda_hat", format = "matrix")

lambda_post_mean <- matrix(NA_real_, N, S)
lambda_post_sd   <- matrix(NA_real_, N, S)
for (s in seq_len(S)) {
  cols <- ((s - 1) * N + 1):(s * N)
  lambda_post_mean[, s] <- colMeans(lambda_hat_draws[, cols])
  lambda_post_sd[, s]   <- apply(lambda_hat_draws[, cols], 2, sd)
}

# Surface samples only (Z_sample == 0) for spatial comparison
surf_idx <- which(as.numeric(samples[["Z_sample"]]) == 0)

# Build plotting data frame — define plot_df here before it is used below
plot_df <- bind_rows(lapply(seq_len(S), function(s) {
  tibble(
    X           = as.numeric(samples[["X"]])[surf_idx],
    Y           = as.numeric(samples[["Y"]])[surf_idx],
    Z_bathy     = as.numeric(samples[["Z_bathy"]])[surf_idx],
    species     = sp_common[s],
    lambda_true = lambda_true_si[surf_idx, s],
    lambda_est  = lambda_post_mean[surf_idx, s],
    lambda_sd   = lambda_post_sd[surf_idx, s]
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

# True vs estimated scatter on log scale
p_scatter <- ggplot(plot_df,
                    aes(x = log(lambda_true + 0.01),
                        y = log(lambda_est  + 0.01),
                        colour = species)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.45, size = 1.4) +
  facet_wrap(~species) +
  scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
  labs(
    title = "True vs posterior mean log-λ (surface samples, Z_sample = 0)",
    x     = "log(true λ)",
    y     = "log(posterior mean λ)"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(OUTPUT_DIR, "lambda_scatter.png"),
       p_scatter, width = 12, height = 4)

# Spatial maps: X vs Y coloured by log posterior mean lambda
sp_colours <- viridis(S, option = "viridis", end = 0.85)

p_maps <- lapply(seq_len(S), function(s) {
  df <- plot_df %>% filter(species == sp_common[s])
  ggplot(df, aes(x = X, y = Y,
                 colour = log(lambda_est + 0.01))) +
    geom_point(size = 2.0, alpha = 0.80) +
    scale_colour_viridis_c(option = "viridis", name = "log(λ̂)") +
    coord_fixed() +
    labs(
      title = sp_common[s],
      x     = "X (km, offshore→nearshore)",
      y     = "Y (km, south→north)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "italic"))
})

ggsave(file.path(OUTPUT_DIR, "lambda_maps.png"),
       wrap_plots(p_maps, nrow = 1),
       width = 15, height = 5)

# R² on log scale per species
for (s in seq_len(S)) {
  df_s <- plot_df %>% filter(species == sp_common[s])
  r2   <- cor(log(df_s$lambda_true + 0.01),
              log(df_s$lambda_est  + 0.01))^2
  cat(sprintf("  %s: R² log-λ = %.3f\n", sp_common[s], r2))
}

# =============================================================================
# 13. Bathymetric preference recovery
# =============================================================================
cat("=== Step 13: Bathymetric preference recovery ===\n")

# plot_df and surf_idx are defined in step 12 above
bathy_recovery_df <- bind_rows(lapply(seq_len(S), function(s) {
  tibble(
    Z_bathy     = as.numeric(samples[["Z_bathy"]])[surf_idx],
    lambda_true = lambda_true_si[surf_idx, s],
    lambda_est  = lambda_post_mean[surf_idx, s],
    species     = sp_common[s]
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

p_bathy_rec <- ggplot(bathy_recovery_df, aes(x = Z_bathy)) +
  geom_point(aes(y = log(lambda_true + 0.01)),
             colour = "grey65", size = 1.0, alpha = 0.4) +
  geom_point(aes(y = log(lambda_est + 0.01), colour = species),
             size = 1.3, alpha = 0.65) +
  geom_smooth(aes(y = log(lambda_true + 0.01)),
              method = "loess", span = 0.4, se = FALSE,
              colour = "grey40", linetype = "dashed", linewidth = 0.8) +
  geom_smooth(aes(y = log(lambda_est + 0.01), colour = species),
              method = "loess", span = 0.4, se = FALSE,
              linewidth = 1.1) +
  facet_wrap(~species) +
  scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
  geom_vline(xintercept = 200, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  annotate("text", x = 215, y = Inf, vjust = 1.4,
           label = "~200 m shelf break", colour = "grey40",
           size = 3.0, hjust = 0, fontface = "italic") +
  labs(
    title    = "Bathymetric preference: true (grey) vs estimated (colour)",
    subtitle = "log(λ) at surface samples; dashed loess = true, solid loess = posterior mean",
    x        = "Bottom depth Z_bathy (m)",
    y        = "log(λ)"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(OUTPUT_DIR, "bathy_preference_recovery.png"),
       p_bathy_rec, width = 12, height = 4)

# =============================================================================
# 14. Save diagnostics and session info
# =============================================================================
cat("=== Step 14: Saving diagnostics ===\n")

write_csv(diag_sum, file.path(OUTPUT_DIR, "diagnostics_summary.csv"))

if (nrow(rhat_bad) > 0)
  write_csv(rhat_bad, file.path(OUTPUT_DIR, "rhat_flagged.csv"))
if (nrow(ess_bad) > 0)
  write_csv(ess_bad,  file.path(OUTPUT_DIR, "ess_flagged.csv"))

sampler_summary <- tibble(
  metric = c("Divergences", "Max treedepth hits",
             "Mean accept stat", "Mean stepsize"),
  value  = c(
    sum(sampler_diag$divergent__),
    sum(sampler_diag$treedepth__ >= MAX_TREEDEPTH),
    mean(sampler_diag$accept_stat__),
    mean(sampler_diag$stepsize__)
  )
)
write_csv(sampler_summary, file.path(OUTPUT_DIR, "sampler_summary.csv"))
print(sampler_summary)

sink(file.path(OUTPUT_DIR, "session_info.txt"))
print(sessionInfo())
sink()

cat("\n=== Done. Outputs written to:", OUTPUT_DIR, "===\n")
cat("Files produced:\n")
for (f in list.files(OUTPUT_DIR)) cat(sprintf("  %s\n", f))

# =============================================================================
# 14. Save diagnostics and session info
# =============================================================================
cat("=== Step 14: Saving diagnostics ===\n")

write_csv(diag_sum, file.path(OUTPUT_DIR, "diagnostics_summary.csv"))

if (nrow(rhat_bad) > 0)
  write_csv(rhat_bad, file.path(OUTPUT_DIR, "rhat_flagged.csv"))
if (nrow(ess_bad) > 0)
  write_csv(ess_bad,  file.path(OUTPUT_DIR, "ess_flagged.csv"))

sampler_summary <- tibble(
  metric = c("Divergences", "Max treedepth hits",
             "Mean accept stat", "Mean stepsize"),
  value  = c(
    sum(sampler_diag$divergent__),
    sum(sampler_diag$treedepth__ >= MAX_TREEDEPTH),
    mean(sampler_diag$accept_stat__),
    mean(sampler_diag$stepsize__)
  )
)
write_csv(sampler_summary, file.path(OUTPUT_DIR, "sampler_summary.csv"))
print(sampler_summary)

sink(file.path(OUTPUT_DIR, "session_info.txt"))
print(sessionInfo())
sink()

cat("\n=== Done. Outputs written to:", OUTPUT_DIR, "===\n")

