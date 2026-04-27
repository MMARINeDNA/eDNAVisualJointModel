# =============================================================================
# 05_check_whale_edna_model_v3.2.R
#
# Diagnostics + posterior predictive checks for the fit produced by
# scripts/04_run_whale_edna_model_v3.2.r.
#
# Inputs:  outputs/whale_edna_output_v3.2/whale_edna_sim_v3.2.rds (truth)
#          outputs/whale_edna_output_v3.2/stan_data.rds
#          outputs/whale_edna_output_v3.2/whale_edna_fit.rds
# Outputs: plots, CSVs, and session info in outputs/whale_edna_output_v3.2/
#
# v3.2: hake-only, qPCR-only. Metabarcoding diagnostics, LOO, and PPCs
# from v3 are removed. kappa is fixed in data, so it is not in the
# parameter recovery plot.
# =============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(loo)
library(patchwork)
library(viridis)
library(dplyr)

set.seed(42)

MAX_TREEDEPTH <- 12

OUTPUT_DIR <- "outputs/whale_edna_output_v3.2"

# -----------------------------------------------------------------------------
# 0. Load fit, simulation truth, and Stan data
# -----------------------------------------------------------------------------
cat("=== Loading fit and context ===\n")
fit       <- readRDS(file.path(OUTPUT_DIR, "whale_edna_fit.rds"))
stan_data <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))
sim       <- readRDS("outputs/whale_edna_output_v3.2/whale_edna_sim_v3.2.rds")

samples         <- sim$design$samples
gp_params       <- sim$truth$gp_params
lambda_true_si  <- sim$truth$lambda_true_si

N         <- stan_data$N
S         <- stan_data$S
sp_common <- sim$meta$sp_common

qpcr_detect_vec <- stan_data$qpcr_detect
qpcr_ct_vec     <- stan_data$qpcr_ct

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
if (nrow(rhat_bad) > 0) print(rhat_bad %>% dplyr::select(variable, Rhat))

sampler_diag <- fit$sampler_diagnostics(format = "df")
cat(sprintf("  Divergences       : %d\n", sum(sampler_diag$divergent__)))
cat(sprintf("  Max treedepth hits: %d\n",
            sum(sampler_diag$treedepth__ >= MAX_TREEDEPTH)))

nuts_df <- bayesplot::nuts_params(fit)
energy_plot <- mcmc_nuts_energy(nuts_df) +
  ggtitle("NUTS Energy")
ggsave(file.path(OUTPUT_DIR, "diag_energy.png"),
       energy_plot, width = 8, height = 4)

# Scalar parameter trace plots. v3.2 sampled scalars: mu_sp, gp_sigma,
# sigma_ct (kappa is fixed in data, gamma/beta_phi are gone).
scalar_pars <- c(
  paste0("mu_sp[",    1:S, "]"),
  paste0("gp_sigma[", 1:S, "]"),
  "sigma_ct"
)

trace_plot <- mcmc_trace(draws_arr, pars = scalar_pars) +
  ggtitle("Trace plots")
ggsave(file.path(OUTPUT_DIR, "diag_trace.png"),
       trace_plot, width = 14, height = 10)

print(fit$summary(scalar_pars))

# -----------------------------------------------------------------------------
# 2. LOO-CV (qPCR only)
# -----------------------------------------------------------------------------
cat("=== LOO-CV ===\n")

loo_qpcr <- loo(fit$draws("log_lik_qpcr", format = "matrix"), cores = 2)
cat("qPCR LOO:\n"); print(loo_qpcr)

# -----------------------------------------------------------------------------
# 3. Parameter recovery
# -----------------------------------------------------------------------------
cat("=== Parameter recovery ===\n")

true_scalar <- tibble(
  variable = c(
    paste0("mu_sp[",    1:S, "]"),
    paste0("gp_sigma[", 1:S, "]"),
    "sigma_ct"
  ),
  true = c(
    sapply(gp_params, `[[`, "mu"),
    sapply(gp_params, `[[`, "sigma"),
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
  labs(title = "Parameter recovery: true vs posterior mean +/- 95% CI",
       x = "True value", y = "Posterior mean") +
  theme_bw(base_size = 11)

ggsave(file.path(OUTPUT_DIR, "param_recovery.png"),
       p_recovery, width = 10, height = 6)

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
       y      = "Posterior mean +/- 95% CI",
       colour = "Species", shape = "Dimension") +
  theme_bw(base_size = 11)

ggsave(file.path(OUTPUT_DIR, "lengthscale_recovery.png"),
       p_ls, width = 9, height = 5)

# -----------------------------------------------------------------------------
# 4. Posterior predictive checks (qPCR only)
# -----------------------------------------------------------------------------
cat("=== Posterior predictive checks ===\n")

pp_qpcr_detect <- fit$draws("pp_qpcr_detect", format = "matrix")
pp_qpcr_ct     <- fit$draws("pp_qpcr_ct",     format = "matrix")

ppc_plots <- list()

ppc_plots[["qpcr_detect"]] <-
  ppc_bars(
    y    = qpcr_detect_vec,
    yrep = pp_qpcr_detect[1:min(200, nrow(pp_qpcr_detect)), ]
  ) +
  ggtitle("PPC: qPCR detection (0/1) - Pacific hake") +
  theme_bw(base_size = 11)

det_idx <- which(qpcr_detect_vec == 1)
if (length(det_idx) >= 5) {
  ppc_plots[["qpcr_ct"]] <-
    ppc_dens_overlay(
      y    = qpcr_ct_vec[det_idx],
      yrep = pp_qpcr_ct[1:min(200, nrow(pp_qpcr_ct)), det_idx, drop = FALSE]
    ) +
    ggtitle("PPC: Ct values (detected only) - Pacific hake") +
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

# lambda_hat is N x S column-major: species s occupies columns
# ((s-1)*N + 1) : (s*N). With S=1 there is just one block.
lambda_hat_draws <- fit$draws("lambda_hat", format = "matrix")

lambda_post_mean <- matrix(NA_real_, N, S)
lambda_post_sd   <- matrix(NA_real_, N, S)
for (s in seq_len(S)) {
  cols <- ((s - 1) * N + 1):(s * N)
  lambda_post_mean[, s] <- colMeans(lambda_hat_draws[, cols])
  lambda_post_sd[, s]   <- apply(lambda_hat_draws[, cols], 2, sd)
}

# All samples are surface (Z_sample = 0) in v3.2
plot_df <- bind_rows(lapply(seq_len(S), function(s) {
  tibble(
    X           = as.numeric(samples[["X"]]),
    Y           = as.numeric(samples[["Y"]]),
    Z_bathy     = as.numeric(samples[["Z_bathy"]]),
    species     = sp_common[s],
    lambda_true = lambda_true_si[, s],
    lambda_est  = lambda_post_mean[, s],
    lambda_sd   = lambda_post_sd[, s]
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
    title = "True vs posterior mean log-lambda (all samples, surface)",
    x     = "log(true lambda)",
    y     = "log(posterior mean lambda)"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(OUTPUT_DIR, "lambda_scatter.png"),
       p_scatter, width = 6, height = 5)

p_maps <- lapply(seq_len(S), function(s) {
  df <- plot_df %>% filter(species == sp_common[s])
  ggplot(df, aes(x = X, y = Y,
                 colour = log(lambda_est + 0.01))) +
    geom_point(size = 2.0, alpha = 0.80) +
    scale_colour_viridis_c(option = "viridis", name = "log(lambda_hat)") +
    coord_fixed() +
    labs(
      title = sp_common[s],
      x     = "X (km, offshore->nearshore)",
      y     = "Y (km, south->north)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "italic"))
})

ggsave(file.path(OUTPUT_DIR, "lambda_maps.png"),
       wrap_plots(p_maps, nrow = 1),
       width = 6, height = 8)

for (s in seq_len(S)) {
  df_s <- plot_df %>% filter(species == sp_common[s])
  r2   <- cor(log(df_s$lambda_true + 0.01),
              log(df_s$lambda_est  + 0.01))^2
  cat(sprintf("  %s: R^2 log-lambda = %.3f\n", sp_common[s], r2))
}

# -----------------------------------------------------------------------------
# 6. Bathymetric preference recovery
# -----------------------------------------------------------------------------
cat("=== Bathymetric preference recovery ===\n")

bathy_recovery_df <- bind_rows(lapply(seq_len(S), function(s) {
  tibble(
    Z_bathy     = as.numeric(samples[["Z_bathy"]]),
    lambda_true = lambda_true_si[, s],
    lambda_est  = lambda_post_mean[, s],
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
    subtitle = "log(lambda) at surface samples; dashed loess = true, solid loess = posterior mean",
    x        = "Bottom depth Z_bathy (m)",
    y        = "log(lambda)"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(OUTPUT_DIR, "bathy_preference_recovery.png"),
       p_bathy_rec, width = 6, height = 5)

# -----------------------------------------------------------------------------
# 7. Posterior decomposition: marginal field structure
#
# Sample log(lambda) = mu_sp + f_s on a regular 3D grid in (X, Y, Z_bathy)
# from each posterior draw, then marginalise to get 1D and 2D posterior
# summaries. This recovers "habitat-preference"-style interpretations
# without fitting structured covariates explicitly. With covariates
# (e.g. log(lambda) = mu + beta_Y * Y_pref(Y) + ... + f_s) the GP basis
# would fight the structured terms because they span the same function
# space; here, the GP does its job (fitting the field) and the
# habitat-preference summaries are read off the posterior afterwards.
# -----------------------------------------------------------------------------
cat("=== Posterior decomposition: marginal field structure ===\n")

# Build the 3D grid spanning the data domain. Limits and HSGP basis size
# are pulled from the saved stan_data so this stays consistent if the
# domain is later changed. We deliberately use a *grid* (not data
# locations) so the marginals are over the full domain, not just where
# stations happen to be.
n_grid_x <- 21
n_grid_y <- 26
n_grid_z <- 36
X_km_max     <- sim$meta$X_km_max
Y_km_max     <- sim$meta$Y_km_max
Z_bathy_max  <- 2 * stan_data$coord_scale[3]   # half-range x 2 = full range

x_grid <- seq(0, X_km_max,    length.out = n_grid_x)
y_grid <- seq(0, Y_km_max,    length.out = n_grid_y)
z_grid <- seq(0, Z_bathy_max, length.out = n_grid_z)

# coord_centre was added to stan_data in a recent format-script revision.
# Older fits (pre-revision) don't have it; in those cases the convention
# was centre = scale (origin at the half-range), so fall back accordingly.
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

# 1-D Laplacian eigenfunctions on [-L, L] (matches stan/whale_edna_hsgp_v3.2.stan)
phi1d_eval <- function(L, j, xn) (1/sqrt(L)) * sin(j * pi * (xn + L) / (2 * L))

phi_x <- sapply(1:m_x, function(j) phi1d_eval(L_hsgp_use[1], j, x_norm_grid))
phi_y <- sapply(1:m_y, function(j) phi1d_eval(L_hsgp_use[2], j, y_norm_grid))
phi_z <- sapply(1:m_z, function(j) phi1d_eval(L_hsgp_use[3], j, z_norm_grid))

# Means over the grid - used to integrate out a dimension when marginalising
phi_x_avg <- colMeans(phi_x)
phi_y_avg <- colMeans(phi_y)
phi_z_avg <- colMeans(phi_z)

# Posterior draws. Thin to ~200 for the marginal computation - cheap enough
# to do the full posterior, but plotting is more legible with a thinned set
# and the credible bands are unchanged.
mu_sp_d_all    <- as.matrix(fit$draws("mu_sp",    format = "matrix"))[, 1]
gp_sigma_d_all <- as.matrix(fit$draws("gp_sigma", format = "matrix"))[, 1]
gp_l_d_all     <- as.matrix(fit$draws("gp_l",     format = "matrix"))   # n x 3 (S=1)
z_beta_d_all   <- as.matrix(fit$draws("z_beta",   format = "matrix"))   # n x M

n_draws_total <- length(mu_sp_d_all)
thin_idx      <- round(seq(1, n_draws_total,
                           length.out = min(200, n_draws_total)))
n_draws_use   <- length(thin_idx)

cat(sprintf("  Grid: %d x %d x %d = %d cells; using %d / %d posterior draws\n",
            n_grid_x, n_grid_y, n_grid_z,
            n_grid_x * n_grid_y * n_grid_z,
            n_draws_use, n_draws_total))

# Storage. f_*_marg holds log(lambda) = mu_sp + f at every grid point /
# marginal axis combination.
f_x_marg  <- matrix(0, nrow = n_draws_use, ncol = n_grid_x)
f_y_marg  <- matrix(0, nrow = n_draws_use, ncol = n_grid_y)
f_z_marg  <- matrix(0, nrow = n_draws_use, ncol = n_grid_z)
f_yz_marg <- array(0, dim = c(n_draws_use, n_grid_y, n_grid_z))

# Constants for spectral weight per dimension (depend on draw via sigma, l)
sqrt_2pi <- sqrt(2 * pi)

for (d_i in seq_along(thin_idx)) {
  d        <- thin_idx[d_i]
  mu_d     <- mu_sp_d_all[d]
  sigma_d  <- gp_sigma_d_all[d]
  lx_d     <- gp_l_d_all[d, 1]
  ly_d     <- gp_l_d_all[d, 2]
  lz_d     <- gp_l_d_all[d, 3]
  zb_d     <- z_beta_d_all[d, ]

  # Per-dim spectral density factors (alpha_dim^2 = sigma^(2/3) per dim)
  alpha_dim <- sigma_d^(1/3)
  freq1 <- (1:m_x) * pi / (2 * L_hsgp_use[1]) / coord_scale_use[1]
  freq2 <- (1:m_y) * pi / (2 * L_hsgp_use[2]) / coord_scale_use[2]
  freq3 <- (1:m_z) * pi / (2 * L_hsgp_use[3]) / coord_scale_use[3]
  spd1  <- alpha_dim^2 * sqrt_2pi * lx_d * exp(-0.5 * (lx_d * freq1)^2)
  spd2  <- alpha_dim^2 * sqrt_2pi * ly_d * exp(-0.5 * (ly_d * freq2)^2)
  spd3  <- alpha_dim^2 * sqrt_2pi * lz_d * exp(-0.5 * (lz_d * freq3)^2)

  # Spectral weights as a 3D array indexed [j1, j2, j3]
  wt_3d <- sqrt(outer(outer(spd1, spd2), spd3))

  # z_beta arrives flattened in Stan's order: col = (j1-1)*m_y*m_z + (j2-1)*m_z + j3
  # i.e. j3 fastest. Reshape and permute to get [j1, j2, j3].
  zb_3d <- aperm(array(zb_d, dim = c(m_z, m_y, m_x)), c(3, 2, 1))
  wz    <- wt_3d * zb_3d

  # 1D marginal f(x): integrate out (y, z)
  #   f(x_i) = sum_{j1,j2,j3} wz[j1,j2,j3] * phi_x[i,j1] * phi_y_avg[j2] * phi_z_avg[j3]
  # Contraction: coeff_x[j1] = sum_{j2,j3} wz[j1,j2,j3] * phi_y_avg[j2] * phi_z_avg[j3]
  coeff_x <- vapply(1:m_x,
                    function(j1) as.numeric(t(phi_y_avg) %*% wz[j1, , ] %*% phi_z_avg),
                    numeric(1))
  f_x_marg[d_i, ] <- mu_d + as.vector(phi_x %*% coeff_x)

  # 1D marginal f(y): integrate out (x, z)
  coeff_y <- vapply(1:m_y,
                    function(j2) as.numeric(t(phi_x_avg) %*% wz[, j2, ] %*% phi_z_avg),
                    numeric(1))
  f_y_marg[d_i, ] <- mu_d + as.vector(phi_y %*% coeff_y)

  # 1D marginal f(z): integrate out (x, y)
  coeff_z <- vapply(1:m_z,
                    function(j3) as.numeric(t(phi_x_avg) %*% wz[, , j3] %*% phi_y_avg),
                    numeric(1))
  f_z_marg[d_i, ] <- mu_d + as.vector(phi_z %*% coeff_z)

  # 2D marginal f(y, z): integrate out x
  # coeff_yz[j2, j3] = sum_{j1} phi_x_avg[j1] * wz[j1, j2, j3]
  coeff_yz <- matrix(0, m_y, m_z)
  for (j1 in 1:m_x) coeff_yz <- coeff_yz + phi_x_avg[j1] * wz[j1, , ]
  f_yz_marg[d_i, , ] <- mu_d + phi_y %*% coeff_yz %*% t(phi_z)
}

# Summarise: posterior median + 50% / 95% credible bands per axis grid
summarise_axis <- function(mat, axis_grid, axis_name) {
  qs <- apply(mat, 2, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
  tibble(
    !!axis_name := axis_grid,
    q025 = qs[1, ], q25 = qs[2, ], q50 = qs[3, ], q75 = qs[4, ], q975 = qs[5, ]
  )
}

f_x_summary <- summarise_axis(f_x_marg, x_grid, "X")
f_y_summary <- summarise_axis(f_y_marg, y_grid, "Y")
f_z_summary <- summarise_axis(f_z_marg, z_grid, "Z_bathy")

# 1D marginal plots
plot_axis_marginal <- function(df, axis_col, axis_lab, title_suffix) {
  ggplot(df, aes(x = .data[[axis_col]], y = q50)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#3B528B", alpha = 0.18) +
    geom_ribbon(aes(ymin = q25,  ymax = q75),  fill = "#3B528B", alpha = 0.32) +
    geom_line(colour = "#3B528B", linewidth = 0.9) +
    labs(
      title    = sprintf("Posterior marginal: log(lambda) vs %s", title_suffix),
      subtitle = sprintf("%s; other dims integrated out over the domain grid. Median + 50%% / 95%% CI.", sp_common[1]),
      x        = axis_lab,
      y        = expression(log(lambda))
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey40"))
}

p_marg_x <- plot_axis_marginal(f_x_summary, "X",       "X (km, easting)",          "Easting")
p_marg_y <- plot_axis_marginal(f_y_summary, "Y",       "Y (km, northing)",         "Latitude")
p_marg_z <- plot_axis_marginal(f_z_summary, "Z_bathy", "Z_bathy (m, bottom depth)", "Bottom depth")

ggsave(file.path(OUTPUT_DIR, "posterior_marginals_1d.png"),
       p_marg_x / p_marg_y / p_marg_z,
       width = 8, height = 10)

# 2D marginal: f(Y, Z_bathy) with X integrated out
yz_post_mean <- apply(f_yz_marg, c(2, 3), mean)
yz_df <- expand.grid(Y = y_grid, Z_bathy = z_grid) %>%
  mutate(log_lambda = as.vector(yz_post_mean))

p_marg_yz <- ggplot(yz_df, aes(x = Y, y = Z_bathy, fill = log_lambda)) +
  geom_raster() +
  scale_y_reverse() +
  scale_fill_viridis_c(option = "viridis", name = expression(log(lambda))) +
  labs(
    title    = "Posterior marginal: log(lambda) over (Y, Z_bathy), X integrated out",
    subtitle = sprintf("%s; posterior mean.", sp_common[1]),
    x        = "Y (km, northing)",
    y        = "Z_bathy (m, bottom depth)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        panel.grid = element_blank())

ggsave(file.path(OUTPUT_DIR, "posterior_marginal_yz.png"),
       p_marg_yz, width = 8, height = 5)

# CSV exports for downstream use
write_csv(f_x_summary, file.path(OUTPUT_DIR, "posterior_marginal_x.csv"))
write_csv(f_y_summary, file.path(OUTPUT_DIR, "posterior_marginal_y.csv"))
write_csv(f_z_summary, file.path(OUTPUT_DIR, "posterior_marginal_z.csv"))

cat(sprintf("  log(lambda) range across 1D marginals: x [%.2f, %.2f], y [%.2f, %.2f], z [%.2f, %.2f]\n",
            min(f_x_summary$q50), max(f_x_summary$q50),
            min(f_y_summary$q50), max(f_y_summary$q50),
            min(f_z_summary$q50), max(f_z_summary$q50)))

# -----------------------------------------------------------------------------
# 8. Save diagnostics and session info
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
