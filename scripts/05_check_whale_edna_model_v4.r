# =============================================================================
# 05_check_whale_edna_model_v4.R
#
# Diagnostics + posterior predictive checks for the fit produced by
# scripts/04_run_whale_edna_model_v4.r.
#
# Inputs:  outputs/whale_edna_sim_v4.rds        (for true parameters / truth)
#          outputs/whale_edna_output_v4/stan_data.rds
#          outputs/whale_edna_output_v4/whale_edna_fit.rds
# Outputs: plots, CSVs, and session info in outputs/whale_edna_output_v4/
# =============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(loo)
library(patchwork)
library(viridis)

set.seed(42)

MAX_TREEDEPTH <- 12    # used for sampler diagnostics reporting

OUTPUT_DIR <- "outputs/whale_edna_output_v4"

# -----------------------------------------------------------------------------
# 0. Load fit, simulation truth, and Stan data
# -----------------------------------------------------------------------------
cat("=== Loading fit and context ===\n")
fit       <- readRDS(file.path(OUTPUT_DIR, "whale_edna_fit.rds"))
stan_data <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))
sim       <- readRDS("outputs/whale_edna_sim_v4.rds")

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

# Energy plot
energy_plot <- mcmc_nuts_energy(sampler_diag) +
  ggtitle("NUTS Energy")
ggsave(file.path(OUTPUT_DIR, "diag_energy.png"),
       energy_plot, width = 8, height = 4)

# Scalar parameter trace plots. alpha_ct and beta_ct are now fixed in
# data and so are not sampled.
scalar_pars <- c(
  paste0("mu_sp[",    1:S, "]"),
  paste0("gp_sigma[", 1:S, "]"),
  "kappa", "sigma_ct",
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
    paste0("gp_sigma[", 1:S, "]"),
    "kappa", "sigma_ct"
  ),
  true = c(
    sapply(gp_params, `[[`, "mu"),
    sapply(gp_params, `[[`, "sigma"),
    sim$truth$qpcr_params$kappa,
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
# 7. Save diagnostics and session info
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
