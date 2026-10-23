# =============================================================================
# fit_visual.r
#
# Fit stan/hsgp_visual.stan to one species of a simulate_visual() dataset and
# summarise convergence + recovery against the simulated truth.
#
#   source("R/fit_visual.r")
#   sim <- simulate_visual(seed = 1)
#   res <- fit_visual(sim, "humpback")
#   res$recovery; res$field
#
# As a script (from the repo root), fits every species and writes
# outputs/visual/<scenario>/:
#   Rscript R/fit_visual.r
# Env: SEED, SCENARIO (output dir name), MX, MY, CHAINS, WARMUP, SAMPLE,
#      ADAPT_DELTA, TREEDEPTH, USE_GP (0/1), USE_BATHY (0/1),
#      SEA_STATE (0/1: simulate segment Beaufort), USE_BFT (0/1; default = SEA_STATE)
# =============================================================================

suppressMessages({
  library(cmdstanr)
  library(posterior)
})
source("R/functions_visual.R")

VISUAL_STAN    <- "stan/hsgp_visual.stan"
VISUAL_INCLUDE <- c("stan/include", file.path(HSGP4EDNA_DIR, "stan", "include"))

compile_visual <- function(...) stan_model_fresh(VISUAL_STAN, VISUAL_INCLUDE, ...)

visual_init <- function(sd) {
  M <- sd$M; G <- sd$use_gp
  rate0 <- max(sum(sd$seg_count) / sum(2 * sd$seg_l * 0.5 * sd$w), 1e-6)  # groups / km^2
  # hazard-rate: start the scale at about w / 2
  ls0 <- if (sd$det_key == 1L) log(sd$w / 2) else sd$log_sigma_prior_mean
  function() {
    init <- list(
    mu_sp = log(rate0 * exp(sd$mu_log_prior_mean)) + rnorm(1, 0, 0.2),
    gp_sigma = array(runif(G, 0.8, 1.2), G),
    gp_l_raw = matrix(runif(G * sd$D1, 0.2, 0.4), G, sd$D1),
    z_beta = matrix(0, G, M), beta_bathy = array(0, sd$K_bathy),
    log_sigma = ls0 + rnorm(1, 0, 0.2),
    beta_size = rnorm(1, 0, 0.005),
    beta_bft = array(rnorm(sd$use_bft_covar, 0, 0.02), sd$use_bft_covar),
    beta_vis = array(rnorm(sd$use_vis_covar, 0, 0.02), sd$use_vis_covar),
    log_b = array(sd$log_b_prior_mean + rnorm(sd$det_key, 0, 0.1), sd$det_key),
    inv_sqrt_phi_nb = array(runif(sd$use_nb, 0.3, 0.7), sd$use_nb),
    mu_s = sd$mu_s_prior_shape / sd$mu_s_prior_rate, phi_s = runif(1, 0.5, 2),
    mu_log_s = sd$mu_log_prior_mean + rnorm(1, 0, 0.1),
    sigma_log_s = sd$sigma_log_prior_shape / sd$sigma_log_prior_rate)
    # Zero-size parameters (use_gp = 0, K_bathy = 0) must be left out: an empty
    # matrix is serialised with the wrong number of dimensions.
    init[lengths(init) > 0]
  }
}

# Fit + summarise. `mod` lets callers reuse one compiled model.
fit_visual <- function(sim, sp, stan_data = NULL, mod = NULL,
                       chains = 4L, warmup = 1000L, sample = 1000L,
                       adapt_delta = 0.9, max_treedepth = 12L, seed = 42L,
                       refresh = 0L, ...) {
  sd  <- if (is.null(stan_data)) format_stan_data_visual(sim, sp, ...) else stan_data
  mod <- if (is.null(mod)) compile_visual() else mod
  t0  <- Sys.time()
  fit <- mod$sample(data = sd, chains = chains, parallel_chains = chains,
                    iter_warmup = warmup, iter_sampling = sample,
                    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                    seed = seed, init = visual_init(sd), refresh = refresh,
                    show_messages = FALSE)
  runtime_min <- as.numeric(Sys.time() - t0, units = "mins")
  c(list(fit = fit, stan_data = sd, runtime_min = runtime_min),
    summarise_visual(fit, sim, sp, sd))
}

summarise_visual <- function(fit, sim, sp, sd) {
  diag <- fit$diagnostic_summary(quiet = TRUE)
  p    <- sim$truth$sp_params[[sp]]
  gp   <- sim$truth$gp_params[[sp]]
  true_ll <- sim$truth$log_lambda[, sp]

  vars <- c("mu_sp", "sigma_c", "beta_size", "mean_group_size", "D_mean", "p_det")
  truth <- c(gp$mu, p$sigma_det,
             if (p$use_size_covar == 1L) p$beta_size_truth else NA,
             p$mean_group_size, mean(exp(true_ll)), NA)
  if (sd$use_gp == 1L) {
    vars  <- c(vars, "gp_sigma[1]", "gp_l[1]", "gp_l[2]")
    truth <- c(truth, gp$sigma, gp$lx, gp$ly)
  }
  if (sd$use_vis_covar == 1L) {
    vars  <- c(vars, "beta_vis[1]")
    truth <- c(truth, if (!is.null(sim$design$segments$vis) && !is.null(p$beta_vis_truth)) p$beta_vis_truth else NA)
  }
  if (sd$use_nb == 1L) {
    vars  <- c(vars, "phi_nb")
    truth <- c(truth, if (is.null(p$nb_size_truth) || !is.finite(p$nb_size_truth)) NA else p$nb_size_truth)
  }
  has_bft <- !is.null(sim$design$segments$bft)   # simulated sea state (truth exists)
  if (sd$use_bft_covar == 1L) {
    vars  <- c(vars, "beta_bft[1]")
    truth <- c(truth, if (has_bft && !is.null(p$beta_bft_truth)) p$beta_bft_truth else NA)
  }
  qs <- summarise_draws(fit$draws(vars), median = stats::median,
                        q025 = ~as.numeric(quantile(.x, 0.025)),
                        q975 = ~as.numeric(quantile(.x, 0.975)),
                        rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  recovery <- data.frame(species = sp, param = vars, truth = truth,
                         median = qs$median, q025 = qs$q025, q975 = qs$q975,
                         rhat = qs$rhat, ess_bulk = qs$ess_bulk)
  recovery$covered <- with(recovery, ifelse(is.na(truth), NA, truth >= q025 & truth <= q975))

  post_ll <- colMeans(fit$draws("log_lambda", format = "draws_matrix"))
  field <- data.frame(species = sp, R2 = cor(post_ll, true_ll)^2,
                      rmse = sqrt(mean((post_ll - true_ll)^2)),
                      bias = mean(post_ll - true_ll))
  # Sea-state signature: slope of the field error on segment Beaufort (log
  # density per Beaufort unit). An unmodelled detection decline reads as low
  # density in rough segments, a negative slope of about beta_bft x
  # dlog(ESW)/dlog(sigma); with the covariate it should be ~0.
  if (has_bft) field$err_bft_slope <- unname(coef(lm(I(post_ll - true_ll) ~ sim$design$segments$bft))[2])
  all_rhat <- fit$summary(c("mu_sp", "gp_sigma", "gp_l_raw", "beta_bathy", "log_sigma",
                            "beta_size", "mu_s", "phi_s", "mu_log_s", "sigma_log_s",
                            if (sd$use_bft_covar == 1L) "beta_bft",
                            if (sd$use_nb == 1L) "inv_sqrt_phi_nb",
                            if (sd$use_vis_covar == 1L) "beta_vis"),
                          "rhat")$rhat
  # Truth-free check that the basis budget is not limiting the fit
  gp_prior <- attr(sd, "gp_prior")
  floor <- if (sd$use_gp == 1L && !is.null(gp_prior))
    cbind(species = sp, floor_check(fit$draws(c("gp_l[1]", "gp_l[2]"), format = "draws_matrix"),
                                    gp_prior)) else NULL
  list(recovery = recovery, field = field, floor = floor, n_det = sd$n, M = sd$M,
       post_log_lambda = post_ll,
       divergences = sum(diag$num_divergent),
       treedepth_hits = sum(diag$num_max_treedepth),
       ebfmi = diag$ebfmi, max_rhat = max(all_rhat, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# Script mode
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  env_int <- function(k, d) as.integer(Sys.getenv(k, d))
  SEED     <- env_int("SEED", "1")
  SCENARIO <- Sys.getenv("SCENARIO", "bathysp_surface")
  MX <- Sys.getenv("MX", ""); MY <- Sys.getenv("MY", "")
  HSGP_M <- if (nzchar(MX) && nzchar(MY)) as.integer(c(MX, MY)) else NULL
  out_dir <- file.path("outputs", "visual", SCENARIO)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  SEA_STATE <- env_int("SEA_STATE", "0") == 1L
  USE_BFT   <- env_int("USE_BFT", as.character(as.integer(SEA_STATE))) == 1L
  sim <- simulate_visual(seed = SEED, sea_state = SEA_STATE)
  saveRDS(sim, file.path(out_dir, "sim.rds"))
  mod <- compile_visual()
  res <- list()
  for (sp in sim$meta$species) {
    cat(sprintf("\n=== %s: %d detections ===\n", sp, nrow(sim$observed[[sp]]$obs)))
    r <- fit_visual(sim, sp, mod = mod,
                    chains = env_int("CHAINS", "4"), warmup = env_int("WARMUP", "1000"),
                    sample = env_int("SAMPLE", "1000"),
                    adapt_delta = as.numeric(Sys.getenv("ADAPT_DELTA", "0.9")),
                    max_treedepth = env_int("TREEDEPTH", "12"),
                    HSGP_M = HSGP_M, use_gp = env_int("USE_GP", "1"),
                    use_bathy = env_int("USE_BATHY", "1") == 1L, use_bft = USE_BFT)
    cat(sprintf("runtime %.1f min | divergences %d | treedepth hits %d | max Rhat %.3f | E-BFMI %s\n",
                r$runtime_min, r$divergences, r$treedepth_hits, r$max_rhat,
                paste(sprintf("%.2f", r$ebfmi), collapse = ", ")))
    print(r$recovery, digits = 3, row.names = FALSE)
    print(r$field, digits = 3, row.names = FALSE)
    if (!is.null(r$floor)) print(r$floor, digits = 3, row.names = FALSE)
    if (any(r$floor$flag)) cat("*** WARNING: length-scale posterior piles up at the basis-budget floor\n")
    if (r$max_rhat > 1.05 || r$divergences > 0)
      cat("*** WARNING: not converged / divergent transitions\n")
    r$fit <- NULL   # CmdStan CSVs live in a temp dir; keep the summaries only
    res[[sp]] <- r
  }
  saveRDS(res, file.path(out_dir, "fit_visual_summary.rds"))
  write.csv(do.call(rbind, lapply(res, `[[`, "recovery")),
            file.path(out_dir, "fit_visual_recovery.csv"), row.names = FALSE)
  cat(sprintf("\nSaved %s/{sim.rds,fit_visual_summary.rds,fit_visual_recovery.csv}\n", out_dir))
}
