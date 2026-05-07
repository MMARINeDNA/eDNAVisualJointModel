# =============================================================================
# 05_check_whale_edna_model_v4.1.R
#
# Diagnostics + posterior predictive checks for the fit produced by
# scripts/04_run_whale_edna_model_v4.1.r.
#
# Inputs:  outputs/whale_edna_output_v4.1/whale_edna_sim_v4.1.rds   (truth)
#          outputs/whale_edna_output_v4.1/stan_data.rds
#          outputs/whale_edna_output_v4.1/whale_edna_fit.rds
# Outputs: plots, CSVs, and session info in outputs/whale_edna_output_v4.1/
#
# Changes from v4 (incorporating v3.2 lessons):
#   * kappa, sigma_ct dropped from scalar trace plots and parameter
#     recovery (fixed in data, not sampled).
#   * MAX_TREEDEPTH = 14 to match the run script.
#   * Posterior decomposition: 1D and 2D marginals of the latent field
#     for each species.
#   * Prior-vs-posterior density overlays for mu_sp, gp_sigma, gp_l
#     across all 3 species.
#   * Per-species posterior predictive plots (already in v4) get
#     additional Z_sample-stratified breakdowns.
# =============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(loo)
library(patchwork)
library(viridis)

set.seed(42)

MAX_TREEDEPTH <- 14    # must match value in 04_run_whale_edna_model_v4.1.r

OUTPUT_DIR <- "outputs/whale_edna_output_v4.1"

# -----------------------------------------------------------------------------
# 0. Load fit, simulation truth, and Stan data
# -----------------------------------------------------------------------------
cat("=== Loading fit and context ===\n")
fit       <- readRDS(file.path(OUTPUT_DIR, "whale_edna_fit.rds"))
stan_data <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))
sim       <- readRDS("outputs/whale_edna_output_v4.1/whale_edna_sim_v4.1.rds")

samples         <- sim$design$samples
gp_params       <- sim$truth$gp_params
lambda_true_si  <- sim$truth$lambda_true_si

N         <- stan_data$N
S         <- stan_data$S
K         <- sim$meta$n_mb_rep
sp_common <- sim$meta$sp_common

qpcr_detect_vec <- stan_data$qpcr_detect
qpcr_ct_vec     <- stan_data$qpcr_ct
mb_reads_long   <- stan_data$mb_reads
N_mb_long       <- stan_data$N_mb_long

# -----------------------------------------------------------------------------
# 1. Diagnostics
# -----------------------------------------------------------------------------
cat("=== Diagnostics ===\n")

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

# Energy plot. bayesplot's mcmc_nuts_* family expects a tidy NUTS data
# frame with columns (Chain, Iteration, Parameter, Value), not the
# CmdStanR sampler_diagnostics(format = "df") shape. Use nuts_params()
# to convert.
nuts_df <- bayesplot::nuts_params(fit)
energy_plot <- mcmc_nuts_energy(nuts_df) +
  ggtitle("NUTS Energy")
ggsave(file.path(OUTPUT_DIR, "diag_energy.png"),
       energy_plot, width = 8, height = 4)

# Scalar parameter trace plots. v4.1: alpha_ct, beta_ct, kappa,
# sigma_ct are all fixed in data and so are not sampled.
scalar_pars <- c(
  paste0("mu_sp[",    1:S, "]"),
  paste0("gp_sigma[", 1:S, "]"),
  paste0("beta0_phi[",  1:S, "]"),
  paste0("gamma0_phi[", 1:S, "]"),
  paste0("gamma1_phi[", 1:S, "]")
)

trace_plot <- mcmc_trace(draws_arr, pars = scalar_pars) +
  ggtitle("Trace plots")
ggsave(file.path(OUTPUT_DIR, "diag_trace.png"),
       trace_plot, width = 14, height = 10)

print(fit$summary(scalar_pars))

# -----------------------------------------------------------------------------
# 2. LOO-CV
# -----------------------------------------------------------------------------
cat("=== LOO-CV ===\n")

loo_qpcr <- loo(fit$draws("log_lik_qpcr", format = "matrix"), cores = 2)
loo_mb   <- loo(fit$draws("log_lik_mb",   format = "matrix"), cores = 2)
cat("qPCR LOO:\n");          print(loo_qpcr)
cat("Metabarcoding LOO:\n"); print(loo_mb)

# -----------------------------------------------------------------------------
# 3. Parameter recovery
# -----------------------------------------------------------------------------
cat("=== Parameter recovery ===\n")

true_scalar <- tibble(
  variable = c(
    paste0("mu_sp[",    1:S, "]"),
    paste0("gp_sigma[", 1:S, "]")
  ),
  true = c(
    sapply(gp_params, `[[`, "mu"),
    sapply(gp_params, `[[`, "sigma")
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

# -----------------------------------------------------------------------------
# 4. Posterior predictive checks
# -----------------------------------------------------------------------------
cat("=== Posterior predictive checks ===\n")

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

# -----------------------------------------------------------------------------
# 5. Spatial lambda: true vs posterior mean
# -----------------------------------------------------------------------------
cat("=== Spatial lambda plots ===\n")

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

# -----------------------------------------------------------------------------
# 6. Bathymetric preference recovery
# -----------------------------------------------------------------------------
cat("=== Bathymetric preference recovery ===\n")

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

# -----------------------------------------------------------------------------
# 7. Posterior decomposition: per-species 1D and 2D marginals of f_s
#
# Sample log(lambda_s) = mu_sp_s + f_s on a 3D grid in (X, Y, Z_bathy)
# from each posterior draw, then marginalise to get per-species 1D and
# 2D field summaries. With S=3 we produce 1D marginals over X, Y,
# Z_bathy stratified by species, plus per-species (Y, Z_bathy)
# heatmaps.
# -----------------------------------------------------------------------------
cat("=== Posterior decomposition: marginal field structure ===\n")

n_grid_x <- 21
n_grid_y <- 26
n_grid_z <- 36

X_km_max     <- sim$meta$X_km_max
Y_km_max     <- sim$meta$Y_km_max
Z_bathy_max  <- 2 * stan_data$coord_scale[3]

x_grid <- seq(0, X_km_max,    length.out = n_grid_x)
y_grid <- seq(0, Y_km_max,    length.out = n_grid_y)
z_grid <- seq(0, Z_bathy_max, length.out = n_grid_z)

coord_centre_use <- if (is.null(stan_data$coord_centre)) {
  as.numeric(stan_data$coord_scale)
} else {
  as.numeric(stan_data$coord_centre)
}
coord_scale_use  <- as.numeric(stan_data$coord_scale)
L_hsgp_use       <- as.numeric(stan_data$L_hsgp)
m_hsgp_use       <- as.integer(stan_data$m_hsgp)
m_x <- m_hsgp_use[1]; m_y <- m_hsgp_use[2]; m_z <- m_hsgp_use[3]
M_total_use      <- m_x * m_y * m_z

x_norm_grid <- (x_grid - coord_centre_use[1]) / coord_scale_use[1]
y_norm_grid <- (y_grid - coord_centre_use[2]) / coord_scale_use[2]
z_norm_grid <- (z_grid - coord_centre_use[3]) / coord_scale_use[3]

stopifnot(all(abs(x_norm_grid) < L_hsgp_use[1]))
stopifnot(all(abs(y_norm_grid) < L_hsgp_use[2]))
stopifnot(all(abs(z_norm_grid) < L_hsgp_use[3]))

# 1-D Laplacian eigenfunctions on [-L, L]
phi1d_eval <- function(L, j, xn) (1/sqrt(L)) * sin(j * pi * (xn + L) / (2 * L))

phi_x <- sapply(1:m_x, function(j) phi1d_eval(L_hsgp_use[1], j, x_norm_grid))
phi_y <- sapply(1:m_y, function(j) phi1d_eval(L_hsgp_use[2], j, y_norm_grid))
phi_z <- sapply(1:m_z, function(j) phi1d_eval(L_hsgp_use[3], j, z_norm_grid))

phi_x_avg <- colMeans(phi_x)
phi_y_avg <- colMeans(phi_y)
phi_z_avg <- colMeans(phi_z)

# Pull posterior draws. mu_sp / gp_sigma / gp_l indexed by species; for
# z_beta we slice out the s-th species's M coefficients per draw using
# Stan's column-major order (z_beta[s,j] -> column (j-1)*S + s).
mu_sp_all    <- as.matrix(fit$draws("mu_sp",    format = "matrix"))   # n_draws x S
gp_sigma_all <- as.matrix(fit$draws("gp_sigma", format = "matrix"))   # n_draws x S
gp_l_all     <- as.matrix(fit$draws("gp_l",     format = "matrix"))   # n_draws x (S*3)
z_beta_all   <- as.matrix(fit$draws("z_beta",   format = "matrix"))   # n_draws x (S*M)

n_draws_total <- nrow(mu_sp_all)
thin_idx      <- round(seq(1, n_draws_total, length.out = min(150, n_draws_total)))
n_draws_use   <- length(thin_idx)

cat(sprintf("  Grid: %d x %d x %d = %d cells; using %d / %d posterior draws; %d species\n",
            n_grid_x, n_grid_y, n_grid_z,
            n_grid_x * n_grid_y * n_grid_z,
            n_draws_use, n_draws_total, S))

# Storage per species
f_x_marg  <- array(0, dim = c(S, n_draws_use, n_grid_x))
f_y_marg  <- array(0, dim = c(S, n_draws_use, n_grid_y))
f_z_marg  <- array(0, dim = c(S, n_draws_use, n_grid_z))
f_yz_marg <- array(0, dim = c(S, n_draws_use, n_grid_y, n_grid_z))

sqrt_2pi <- sqrt(2 * pi)

# Stan column index for z_beta[s, j]: (j-1)*S + s.
z_beta_cols <- function(s) ( (1:M_total_use - 1L) * S + s )
gp_l_cols   <- function(s) c(s, S + s, 2L * S + s)   # gp_l[s, 1..3]

for (s in seq_len(S)) {
  zb_cols  <- z_beta_cols(s)
  gpl_cols <- gp_l_cols(s)

  for (d_i in seq_along(thin_idx)) {
    d        <- thin_idx[d_i]
    mu_d     <- mu_sp_all[d, s]
    sigma_d  <- gp_sigma_all[d, s]
    lx_d     <- gp_l_all[d, gpl_cols[1]]
    ly_d     <- gp_l_all[d, gpl_cols[2]]
    lz_d     <- gp_l_all[d, gpl_cols[3]]
    zb_d     <- z_beta_all[d, zb_cols]

    alpha_dim <- sigma_d^(1/3)
    freq1 <- (1:m_x) * pi / (2 * L_hsgp_use[1]) / coord_scale_use[1]
    freq2 <- (1:m_y) * pi / (2 * L_hsgp_use[2]) / coord_scale_use[2]
    freq3 <- (1:m_z) * pi / (2 * L_hsgp_use[3]) / coord_scale_use[3]
    spd1  <- alpha_dim^2 * sqrt_2pi * lx_d * exp(-0.5 * (lx_d * freq1)^2)
    spd2  <- alpha_dim^2 * sqrt_2pi * ly_d * exp(-0.5 * (ly_d * freq2)^2)
    spd3  <- alpha_dim^2 * sqrt_2pi * lz_d * exp(-0.5 * (lz_d * freq3)^2)
    wt_3d <- sqrt(outer(outer(spd1, spd2), spd3))

    zb_3d <- aperm(array(zb_d, dim = c(m_z, m_y, m_x)), c(3, 2, 1))
    wz    <- wt_3d * zb_3d

    coeff_x <- vapply(1:m_x,
                      function(j1) as.numeric(t(phi_y_avg) %*% wz[j1, , ] %*% phi_z_avg),
                      numeric(1))
    f_x_marg[s, d_i, ] <- mu_d + as.vector(phi_x %*% coeff_x)

    coeff_y <- vapply(1:m_y,
                      function(j2) as.numeric(t(phi_x_avg) %*% wz[, j2, ] %*% phi_z_avg),
                      numeric(1))
    f_y_marg[s, d_i, ] <- mu_d + as.vector(phi_y %*% coeff_y)

    coeff_z <- vapply(1:m_z,
                      function(j3) as.numeric(t(phi_x_avg) %*% wz[, , j3] %*% phi_y_avg),
                      numeric(1))
    f_z_marg[s, d_i, ] <- mu_d + as.vector(phi_z %*% coeff_z)

    coeff_yz <- matrix(0, m_y, m_z)
    for (j1 in 1:m_x) coeff_yz <- coeff_yz + phi_x_avg[j1] * wz[j1, , ]
    f_yz_marg[s, d_i, , ] <- mu_d + phi_y %*% coeff_yz %*% t(phi_z)
  }
}

# 1D marginal summaries — long form for ggplot
summarise_axis_multispecies <- function(arr, axis_grid, axis_name) {
  bind_rows(lapply(seq_len(S), function(s) {
    qs <- apply(arr[s, , , drop = FALSE], 3, quantile,
                probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
    tibble(
      species = sp_common[s],
      !!axis_name := axis_grid,
      q025 = qs[1, ], q25 = qs[2, ], q50 = qs[3, ],
      q75  = qs[4, ], q975 = qs[5, ]
    )
  }))
}

f_x_summary <- summarise_axis_multispecies(f_x_marg, x_grid, "X")       %>%
                  mutate(species = factor(species, levels = sp_common))
f_y_summary <- summarise_axis_multispecies(f_y_marg, y_grid, "Y")       %>%
                  mutate(species = factor(species, levels = sp_common))
f_z_summary <- summarise_axis_multispecies(f_z_marg, z_grid, "Z_bathy") %>%
                  mutate(species = factor(species, levels = sp_common))

plot_axis_marginal_multi <- function(df, axis_col, axis_lab) {
  ggplot(df, aes(x = .data[[axis_col]], y = q50,
                 colour = species, fill = species)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.18, colour = NA) +
    geom_ribbon(aes(ymin = q25,  ymax = q75),  alpha = 0.32, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~species, scales = "free_y") +
    scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
    scale_fill_viridis_d  (option = "D", end = 0.85, guide = "none") +
    labs(x = axis_lab, y = expression(log(lambda))) +
    theme_bw(base_size = 11)
}

p_marg_x <- plot_axis_marginal_multi(f_x_summary, "X",       "X (km, easting)")
p_marg_y <- plot_axis_marginal_multi(f_y_summary, "Y",       "Y (km, northing)")
p_marg_z <- plot_axis_marginal_multi(f_z_summary, "Z_bathy", "Z_bathy (m, bottom depth)")

ggsave(file.path(OUTPUT_DIR, "posterior_marginals_1d.png"),
       p_marg_x / p_marg_y / p_marg_z +
         plot_annotation(title = "Per-species posterior marginals: log(lambda)",
                         subtitle = "Median + 50% / 95% CI; other dims integrated out"),
       width = 12, height = 12)

# 2D marginals: posterior mean f_yz per species
yz_post_mean <- array(0, dim = c(S, n_grid_y, n_grid_z))
for (s in seq_len(S)) {
  yz_post_mean[s, , ] <- apply(f_yz_marg[s, , , , drop = FALSE], c(3, 4), mean)
}
yz_df <- bind_rows(lapply(seq_len(S), function(s) {
  expand.grid(Y = y_grid, Z_bathy = z_grid) %>%
    mutate(log_lambda = as.vector(yz_post_mean[s, , ]),
           species = sp_common[s])
})) %>% mutate(species = factor(species, levels = sp_common))

p_marg_yz <- ggplot(yz_df, aes(x = Y, y = Z_bathy, fill = log_lambda)) +
  geom_raster() +
  facet_wrap(~species, scales = "free") +
  scale_y_reverse() +
  scale_fill_viridis_c(option = "viridis", name = expression(log(lambda))) +
  labs(title = "Posterior marginal: log(lambda) over (Y, Z_bathy), X integrated out",
       x = "Y (km, northing)", y = "Z_bathy (m, bottom depth)") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

ggsave(file.path(OUTPUT_DIR, "posterior_marginal_yz.png"),
       p_marg_yz, width = 14, height = 5)

write_csv(f_x_summary, file.path(OUTPUT_DIR, "posterior_marginal_x.csv"))
write_csv(f_y_summary, file.path(OUTPUT_DIR, "posterior_marginal_y.csv"))
write_csv(f_z_summary, file.path(OUTPUT_DIR, "posterior_marginal_z.csv"))

# -----------------------------------------------------------------------------
# 8. Prior vs posterior densities, per species
#
# For each species and each scalar GP parameter (mu_sp, gp_sigma, lx,
# ly, lz), overlay the prior density on a kernel-density estimate of
# the posterior, with a vertical dashed line at the simulated truth.
# All gp_l priors are identical across species in v4.1; what differs is
# the truth (the species' kernel scales) and the posterior.
# -----------------------------------------------------------------------------
cat("=== Prior vs posterior density plots ===\n")

dnorm_trunc0 <- function(x, mu, sig) {
  ifelse(x >= 0, dnorm(x, mu, sig) / (1 - pnorm(0, mu, sig)), 0)
}

# Per-species spec list. Each row -> one panel.
prior_post_specs <- bind_rows(lapply(seq_len(S), function(s) {
  tibble::tribble(
    ~var,                    ~label,                        ~truth,                              ~prior_kind, ~p1, ~p2,
    sprintf("mu_sp[%d]", s),    "mu_sp",                  gp_params[[s]]$mu,                    "norm", stan_data$prior_mu_sp_mu,    stan_data$prior_mu_sp_sig,
    sprintf("gp_sigma[%d]", s), "gp_sigma",               gp_params[[s]]$sigma,                 "gamma",stan_data$prior_gp_sigma_shape, stan_data$prior_gp_sigma_rate,
    sprintf("gp_l[%d,1]", s),   "lx (km)",                gp_params[[s]]$lx,                    "trunc",stan_data$prior_gp_lx_mu,    stan_data$prior_gp_lx_sig,
    sprintf("gp_l[%d,2]", s),   "ly (km)",                gp_params[[s]]$ly,                    "trunc",stan_data$prior_gp_ly_mu,    stan_data$prior_gp_ly_sig,
    sprintf("gp_l[%d,3]", s),   "lz_bathy (m)",           gp_params[[s]]$lz,                    "trunc",stan_data$prior_gp_lz_mu,    stan_data$prior_gp_lz_sig
  ) %>% mutate(species = sp_common[s])
})) %>% mutate(species = factor(species, levels = sp_common))

prior_density <- function(kind, p1, p2, x) {
  switch(kind,
    norm  = dnorm(x, p1, p2),
    gamma = dgamma(x, p1, p2),
    trunc = dnorm_trunc0(x, p1, p2),
    stop("unknown prior_kind"))
}
prior_xrange <- function(kind, p1, p2) {
  switch(kind,
    norm  = p1 + c(-3, 3) * p2,
    gamma = c(0, qgamma(0.995, p1, p2)),
    trunc = c(0, p1 + 3 * p2),
    stop("unknown prior_kind"))
}

prior_post_plots <- lapply(seq_len(nrow(prior_post_specs)), function(k) {
  spec  <- prior_post_specs[k, ]
  rng   <- prior_xrange(spec$prior_kind, spec$p1, spec$p2)
  draws_v <- as.numeric(fit$draws(spec$var, format = "matrix"))

  post_dens <- density(draws_v,
                       from = max(min(rng), min(draws_v) - 0.05 * diff(rng)),
                       to   = min(max(rng), max(draws_v) + 0.05 * diff(rng)),
                       n    = 512)
  post_df  <- tibble(x = post_dens$x, y = post_dens$y, dist = "posterior")
  x_grid_p <- seq(min(rng), max(rng), length.out = 512)
  prior_df <- tibble(x = x_grid_p,
                     y = prior_density(spec$prior_kind, spec$p1, spec$p2, x_grid_p),
                     dist = "prior")

  ggplot(bind_rows(post_df, prior_df),
         aes(x = x, y = y, colour = dist, fill = dist)) +
    geom_area(alpha = 0.30, position = "identity", colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = spec$truth, linetype = "dashed",
               colour = "red", linewidth = 0.5) +
    scale_colour_manual(values = c(posterior = "#3B528B", prior = "grey40")) +
    scale_fill_manual  (values = c(posterior = "#3B528B", prior = "grey60")) +
    labs(title = sprintf("%s -- %s", spec$species, spec$label),
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(plot.title    = element_text(face = "bold", size = 9),
          legend.position = "none",
          axis.text     = element_text(size = 7))
})

# 5 columns (parameters) x S rows (species)
prior_post_combined <- wrap_plots(prior_post_plots, ncol = 5) +
  plot_annotation(
    title    = "Prior vs posterior density (per-species)",
    subtitle = "Vertical dashed line = simulated truth.  Columns: mu_sp, gp_sigma, lx, ly, lz.",
    theme    = theme(plot.title = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 10, colour = "grey40"))
  )

ggsave(file.path(OUTPUT_DIR, "prior_vs_posterior.png"),
       prior_post_combined, width = 16, height = 3 * S + 1)

# -----------------------------------------------------------------------------
# 9. Save diagnostics and session info
# -----------------------------------------------------------------------------
cat("=== Saving diagnostics ===\n")

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
