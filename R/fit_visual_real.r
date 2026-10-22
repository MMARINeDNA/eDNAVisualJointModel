# =============================================================================
# fit_visual_real.r
#
# Fit stan/hsgp_visual.stan to the real 35N line-transect data
# (data/processed/muri_visual.rds) for one species, with the detection,
# sea-state, visibility and count settings stored in the object
# (truth$sp_params: det_key, use_size_covar, use_bft, use_vis, use_nb), and
# save a summary: diagnostics, the length-scale floor check, the main
# parameters and the posterior mean field at the segments.
#
#   SPECIES=pwsd Rscript R/fit_visual_real.r
# Env: SPECIES (humpback | pwsd), M_MAX (basis budget; GP_BASIS_BUDGET = 800),
#      CHAINS, WARMUP, SAMPLE (4 x 1000/1000), ADAPT_DELTA (0.9), TREEDEPTH (12),
#      TAG (output name; default "m<M_MAX>"), OUT (outputs/visual/real)
# Out: <OUT>/<species>_<TAG>.rds
# =============================================================================

source("R/fit_visual.r")

env_int <- function(k, d) as.integer(Sys.getenv(k, d))
sp    <- Sys.getenv("SPECIES", "humpback")
M_MAX <- env_int("M_MAX", as.character(GP_BASIS_BUDGET))
TAG   <- Sys.getenv("TAG", sprintf("m%d", M_MAX))
OUT   <- Sys.getenv("OUT", "outputs/visual/real")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

vis <- readRDS("data/processed/muri_visual.rds")
sd  <- format_stan_data_visual(vis, sp, M_max = M_MAX)
cat(sprintf("%s: n = %d, n_seg = %d, w = %.2f, det_key = %d, use_bft = %d, use_vis = %d, use_nb = %d, cells = %d, M = %d\n",
            sp, sd$n, sd$n_seg, sd$w, sd$det_key, sd$use_bft_covar, sd$use_vis_covar, sd$use_nb,
            sd$n_bft, sd$M))

mod <- compile_visual()
t0  <- Sys.time()
fit <- mod$sample(data = sd, chains = env_int("CHAINS", "4"), parallel_chains = env_int("CHAINS", "4"),
                  iter_warmup = env_int("WARMUP", "1000"), iter_sampling = env_int("SAMPLE", "1000"),
                  adapt_delta = as.numeric(Sys.getenv("ADAPT_DELTA", "0.9")),
                  max_treedepth = env_int("TREEDEPTH", "12"), seed = 42L, init = visual_init(sd),
                  refresh = 0, show_messages = FALSE)
runtime <- as.numeric(Sys.time() - t0, units = "mins")

vars <- c("mu_sp", "gp_sigma[1]", "gp_l[1]", "gp_l[2]", "D_mean", "p_det", "sigma_c",
          "mean_group_size", if (sd$det_key == 1L) "b_shape", if (sd$use_size_covar == 1L) "beta_size",
          if (sd$use_bft_covar == 1L) "beta_bft[1]", if (sd$use_vis_covar == 1L) "beta_vis[1]",
          if (sd$use_nb == 1L) "phi_nb")
summ <- fit$summary(vars, median = stats::median, q05 = ~quantile(.x, 0.05), q95 = ~quantile(.x, 0.95),
                    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
diag <- fit$diagnostic_summary(quiet = TRUE)
max_rhat <- max(fit$summary(c("mu_sp", "gp_sigma", "gp_l_raw", "beta_bathy", "log_sigma", "mu_s",
                              "mu_log_s", if (sd$use_bft_covar == 1L) "beta_bft",
                              if (sd$use_vis_covar == 1L) "beta_vis",
                              if (sd$use_nb == 1L) "inv_sqrt_phi_nb"), "rhat")$rhat, na.rm = TRUE)
floor <- floor_check(fit$draws(c("gp_l[1]", "gp_l[2]"), format = "draws_matrix"), attr(sd, "gp_prior"))
res <- list(species = sp, tag = TAG, M = sd$M, M_max = M_MAX, n = sd$n, n_cells = sd$n_bft,
            settings = sd[c("det_key", "use_size_covar", "use_bft_covar", "use_vis_covar", "use_nb")],
            runtime_min = runtime, divergences = sum(diag$num_divergent),
            treedepth_hits = sum(diag$num_max_treedepth), ebfmi = diag$ebfmi, max_rhat = max_rhat,
            summary = summ, floor = floor,
            post_log_lambda = colMeans(fit$draws("log_lambda", format = "draws_matrix")))
f <- file.path(OUT, sprintf("%s_%s.rds", sp, TAG))
saveRDS(res, f)
cat(sprintf("runtime %.1f min | divergences %d | treedepth hits %d | max Rhat %.3f | E-BFMI %s\n",
            runtime, res$divergences, res$treedepth_hits, max_rhat, paste(sprintf("%.2f", diag$ebfmi), collapse = ", ")))
print(as.data.frame(summ), digits = 3, row.names = FALSE)
print(floor, digits = 3, row.names = FALSE)
cat("Saved", f, "\n")
