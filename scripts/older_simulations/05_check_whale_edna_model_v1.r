# =============================================================================
# 05_check_whale_edna_model_v1.R
#
# Diagnostics + posterior predictive checks for the fit produced by
# scripts/04_run_whale_edna_model_v1.r.
#
# Inputs:  outputs/whale_edna_sim_v1.rds
#          outputs/whale_edna_output_v1/stan_data.rds
#          outputs/whale_edna_output_v1/whale_edna_fit.rds
# Outputs: plots, CSVs, and session info in outputs/whale_edna_output_v1/
# =============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(loo)
library(patchwork)
library(viridis)

set.seed(42)

MAX_TREEDEPTH <- 12
OUTPUT_DIR    <- "outputs/whale_edna_output_v1"

cat("=== Loading fit and context ===\n")
fit       <- readRDS(file.path(OUTPUT_DIR, "whale_edna_fit.rds"))
stan_data <- readRDS(file.path(OUTPUT_DIR, "stan_data.rds"))
sim       <- readRDS("outputs/whale_edna_sim_v1.rds")

samples       <- sim$design$samples
lambda_true   <- sim$truth$lambda_true
sp_common     <- sim$meta$sp_common

N <- stan_data$N
S <- stan_data$S
K <- stan_data$K

n_detect     <- stan_data$n_detect
any_detect   <- stan_data$any_detect
mean_ct_obs  <- stan_data$mean_ct_obs
mb_reads_arr <- stan_data$mb_reads   # N × S × K
mb_total_arr <- stan_data$mb_total   # N × K

sp_colours <- viridis(S, option = "viridis", end = 0.85)

# -----------------------------------------------------------------------------
# 1. Diagnostics
# -----------------------------------------------------------------------------
cat("=== Diagnostics ===\n")
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

# bayesplot's mcmc_nuts_* family expects a tidy NUTS data frame (Chain,
# Iteration, Parameter, Value), not the wide CmdStanR shape we use for
# the divergence / treedepth counters above. Use nuts_params() to get
# the right format directly from the fit.
nuts_df <- bayesplot::nuts_params(fit)
ggsave(file.path(OUTPUT_DIR, "diag_energy.png"),
       mcmc_nuts_energy(nuts_df) +
         ggtitle("NUTS Energy — Bayesian Fraction of Missing Information"),
       width = 8, height = 4)

scalar_pars <- c(
  paste0("gp_sigma[", 1:S, "]"),
  paste0("mu_sp[",    1:S, "]"),
  "kappa",
  "p_zi",
  paste0("phi_bb[",   1:S, "]")
)
ggsave(file.path(OUTPUT_DIR, "diag_trace.png"),
       mcmc_trace(draws_arr, pars = scalar_pars) +
         ggtitle("Trace plots — scalar parameters"),
       width = 12, height = 10)

# -----------------------------------------------------------------------------
# 2. LOO-CV
# -----------------------------------------------------------------------------
cat("=== LOO-CV ===\n")
log_lik_qpcr <- fit$draws("log_lik_qpcr", format = "matrix")
log_lik_mb   <- fit$draws("log_lik_mb",   format = "matrix")
loo_qpcr <- loo(log_lik_qpcr, cores = 2)
loo_mb   <- loo(log_lik_mb,   cores = 2)
cat("  LOO — qPCR:\n");          print(loo_qpcr)
cat("  LOO — Metabarcoding:\n"); print(loo_mb)

png(file.path(OUTPUT_DIR, "loo_pareto_k.png"), width = 900, height = 400)
par(mfrow = c(1, 2))
plot(loo_qpcr, main = "LOO Pareto-k: qPCR (hake)")
plot(loo_mb,   main = "LOO Pareto-k: Metabarcoding")
dev.off()

# -----------------------------------------------------------------------------
# 3. Parameter recovery
# -----------------------------------------------------------------------------
cat("=== Parameter recovery ===\n")

scalar_summary <- diag_summary %>%
  filter(variable %in% scalar_pars) %>%
  mutate(
    param  = str_extract(variable, "^[a-z_]+"),
    sp_idx = as.integer(str_extract(variable, "[0-9]+"))
  )

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
  mutate(sp_name = ifelse(is.na(sp_idx), "hake", sp_common[sp_idx]))

p_recovery <- ggplot(recovery_df,
                     aes(x = true, y = mean,
                         ymin = `2.5%`, ymax = `97.5%`,
                         colour = sp_name, shape = sp_name)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(width = 0, alpha = 0.6) +
  geom_point(size = 3) +
  facet_wrap(~param, scales = "free") +
  scale_colour_viridis_d(option = "D", end = 0.85) +
  labs(title  = "Parameter recovery: true vs posterior mean (95% CI)",
       x = "True value", y = "Posterior mean ± 95% CI",
       colour = "Species", shape = "Species") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTPUT_DIR, "param_recovery.png"), p_recovery, width = 12, height = 8)

ls_pars <- paste0("gp_l[", rep(1:S, each = 3), ",", rep(1:3, S), "]")
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
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(width = 0, alpha = 0.6) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "C", end = 0.85) +
  labs(title = "Length-scale recovery (lx, ly in km; lz_bathy in m)",
       x = "True length-scale", y = "Posterior mean ± 95% CI",
       colour = "Species", shape = "Dimension") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTPUT_DIR, "lengthscale_recovery.png"), p_ls, width = 8, height = 5)

# -----------------------------------------------------------------------------
# 4. Posterior predictive checks
# -----------------------------------------------------------------------------
cat("=== Posterior predictive checks ===\n")
pp_n_detect <- fit$draws("pp_n_detect", format = "matrix")
pp_mean_ct  <- fit$draws("pp_mean_ct",  format = "matrix")
pp_mb_reads <- fit$draws("pp_mb_reads", format = "matrix")

# Stan stores [N, S] column-major: species s column block = ((s-1)*N + 1):(s*N)
reshape_ns <- function(mat, N, S) {
  lapply(1:S, function(s) mat[, ((s-1)*N + 1):(s*N)])
}
# mb_reads is [N, S, K]: Stan column-major over all three dims
reshape_nsk <- function(mat, N, S, K) {
  lapply(1:S, function(s) {
    lapply(1:K, function(k) {
      col_start <- ((k-1)*S*N) + ((s-1)*N) + 1
      mat[, col_start:(col_start + N - 1)]
    })
  })
}

pp_det_list <- list(pp_n_detect)
pp_ct_list  <- list(pp_mean_ct)
pp_mb_list  <- reshape_nsk(pp_mb_reads, N, S, K)

ppc_plots <- list()

# qPCR (hake only)
ppc_plots[["det_hake"]] <-
  ppc_bars(n_detect, pp_det_list[[1]][1:200, ]) +
  ggtitle("PPC: qPCR n_detect — Pacific hake") +
  theme_bw(base_size = 11)

idx_det <- which(any_detect == 1)
if (length(idx_det) > 5) {
  ppc_plots[["ct_hake"]] <-
    ppc_dens_overlay(mean_ct_obs[idx_det],
                     pp_ct_list[[1]][1:200, idx_det]) +
    ggtitle("PPC: mean Ct (detected only) — Pacific hake") +
    theme_bw(base_size = 11)
}

for (s in seq_len(S)) {
  y_mb <- mb_reads_arr[, s, 1]
  ppc_plots[[paste0("mb_", s)]] <-
    ppc_dens_overlay(y_mb, pp_mb_list[[s]][[1]][1:100, ]) +
    ggtitle(sprintf("PPC: MB reads rep 1 — %s", sp_common[s])) +
    theme_bw(base_size = 11)
}

ggsave(file.path(OUTPUT_DIR, "ppc_all.png"),
       wrap_plots(ppc_plots, ncol = 2),
       width = 12, height = 4 * ceiling(length(ppc_plots) / 2))

# -----------------------------------------------------------------------------
# 5. Spatial lambda: true vs estimated
# -----------------------------------------------------------------------------
cat("=== Spatial lambda plots ===\n")
lambda_hat_draws <- fit$draws("lambda_hat", format = "matrix")
lambda_hat_list  <- reshape_ns(lambda_hat_draws, N, S)

lambda_post_mean <- matrix(NA, N, S)
lambda_post_sd   <- matrix(NA, N, S)
for (s in seq_len(S)) {
  lambda_post_mean[, s] <- colMeans(lambda_hat_list[[s]])
  lambda_post_sd[, s]   <- apply(lambda_hat_list[[s]], 2, sd)
}

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
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.4, size = 1.2) +
  facet_wrap(~species) +
  scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
  labs(title = "True vs posterior mean log-lambda (surface samples)",
       x = "log(true λ)", y = "log(posterior mean λ)") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTPUT_DIR, "lambda_scatter.png"), p_scatter, width = 12, height = 4)

# -----------------------------------------------------------------------------
# 6. Bathymetric preference recovery
# -----------------------------------------------------------------------------
cat("=== Bathymetric preference recovery ===\n")
bathy_recovery_df <- plot_df   # reuse same table

p_bathy_rec <- ggplot(bathy_recovery_df, aes(x = Z_bathy)) +
  geom_point(aes(y = log(lambda_true + 0.01)),
             colour = "grey60", size = 1.0, alpha = 0.4) +
  geom_point(aes(y = log(lambda_est + 0.01), colour = species),
             size = 1.2, alpha = 0.6) +
  facet_wrap(~species) +
  scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
  labs(title    = "Bathymetric preference: true (grey) vs estimated (colour)",
       subtitle = "log(λ) vs Z_bathy at surface samples; grey = true, colour = posterior mean",
       x = "Bottom depth Z_bathy (m)", y = "log(λ)") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTPUT_DIR, "bathy_preference_recovery.png"),
       p_bathy_rec, width = 12, height = 4)

# True bathymetric mean function
zbathy_seq <- seq(10, 2500, by = 20)

zbathy_pref_fn <- function(Z_bathy, mu_z, sd_z, amp) {
  raw <- dnorm(Z_bathy, mu_z, sd_z)
  raw <- raw / max(raw)
  amp * (raw - 0.5)
}

bathy_fn_df <- bind_rows(lapply(seq_len(S), function(s) {
  p <- sim$truth$gp_params[[s]]
  z_mu_true <- zbathy_pref_fn(zbathy_seq,
                              p$zbathy_pref_mu,
                              p$zbathy_pref_sd,
                              p$zbathy_pref_amp)
  data.frame(
    Z_bathy    = zbathy_seq,
    log_lambda = p$mu + z_mu_true,
    type       = "True mean function",
    species    = sp_common[s]
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

p_bathy_fn <- ggplot() +
  geom_line(data      = bathy_fn_df,
            aes(x = Z_bathy, y = log_lambda, colour = species),
            linewidth = 1.2) +
  geom_vline(xintercept = 200, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  annotate("text", x = 210, y = max(bathy_fn_df$log_lambda) * 0.98,
           label = "shelf break\n(~200 m)", colour = "grey40",
           size = 3.0, hjust = 0, fontface = "italic") +
  scale_colour_manual(values = sp_colours, name = NULL, labels = sp_common) +
  scale_x_continuous(
    breaks = c(0, 200, 500, 1000, 1500, 2000, 2500),
    expand = c(0.02, 0)
  ) +
  labs(title    = "True bathymetric mean function by species",
       subtitle = "log(λ) = mu_s + zbathy_pref(Z_bathy)",
       x = "Bottom depth Z_bathy (m)", y = "log(λ)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        legend.text     = element_text(face = "italic", size = 9))
ggsave(file.path(OUTPUT_DIR, "bathy_mean_function.png"),
       p_bathy_fn, width = 9, height = 5)

# -----------------------------------------------------------------------------
# 7. Save diagnostics summary
# -----------------------------------------------------------------------------
cat("=== Saving diagnostics ===\n")
write_csv(diag_summary, file.path(OUTPUT_DIR, "diagnostics_summary.csv"))
if (nrow(rhat_bad) > 0) write_csv(rhat_bad, file.path(OUTPUT_DIR, "rhat_flagged.csv"))
if (nrow(ess_bad)  > 0) write_csv(ess_bad,  file.path(OUTPUT_DIR, "ess_flagged.csv"))

sampler_summary <- data.frame(
  metric = c("Divergences", "Max treedepth hits",
             "Mean accept stat", "Mean stepsize"),
  value  = c(
    n_div, n_max,
    mean(sampler_diag$accept_stat__),
    mean(sampler_diag$stepsize__)
  )
)
write_csv(sampler_summary, file.path(OUTPUT_DIR, "sampler_summary.csv"))
print(sampler_summary)

sink(file.path(OUTPUT_DIR, "session_info.txt"))
print(sessionInfo())
sink()

cat("\n=== All done. Outputs written to:", OUTPUT_DIR, "===\n")
cat("Files produced:\n")
for (f in list.files(OUTPUT_DIR)) cat(sprintf("  %s\n", f))
