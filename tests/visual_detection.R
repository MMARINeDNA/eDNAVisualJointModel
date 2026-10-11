# =============================================================================
# tests/visual_detection.R
#
# Checks the detection-function options of stan/hsgp_visual.stan on simulated
# data (non-spatial: use_gp = 0, so it is quick):
#
#   humpback-like   hazard-rate key + group-size covariate
#   pwsd-like       half-normal key + group size + two segment covariates
#
# For each: (1) the true detection parameters lie inside the 99% posterior
# interval, and (2) Stan's per-covariate-row ESW agrees with R's integrate()
# (the hazard-rate ESW is a midpoint rule; the group-size sum is binned).
#
#   Rscript tests/visual_detection.R          # from the repo root
#
# NOTE: cmdstanr does not notice edits to stan/include/*.stan; the model is
# compiled with force_recompile = TRUE here for that reason.
# =============================================================================

suppressMessages({ library(cmdstanr); library(posterior) })
source("R/fit_visual.r")

mod <- cmdstan_model(VISUAL_STAN, include_paths = VISUAL_INCLUDE, force_recompile = TRUE)

sp <- lt_species_params()
sp$humpback <- modifyList(sp$humpback, list(det_key = 1L, sigma_det = 2.5, hr_b_truth = 3,
                                            use_size_covar = 1L, beta_size_truth = 0.1))
sp$pwsd <- modifyList(sp$pwsd, list(det_key = 0L, sigma_det = 1.5,
                                    det_covars = c("avgvis", "avgbft"),
                                    beta_det_truth = c(0.3, -0.3), beta_size_truth = 0.005))
gp <- default_gp_params(); gp$pwsd$mu <- log(0.3); gp$humpback$mu <- log(0.02)
sim <- simulate_visual(seed = 11, n_transects = 25, seg_length = 10, w = 6.5,
                       gp_params = gp, sp_params = sp)

esw_rel_diff <- function(fit, sd) {
  vars <- c("log_sigma", "beta_size", if (sd$det_key == 1L) "log_b",
            if (sd$K_det > 0) "beta_det", "pk", "esw_zc")
  d <- as.data.frame(fit$draws(vars, format = "draws_df"))
  out <- numeric(0)
  for (i in 1:5) {
    pk <- as.numeric(d[i, grep("^pk\\[", names(d))]); K <- length(pk)
    bd <- as.numeric(d[i, grep("^beta_det\\[", names(d))])
    b  <- if (sd$det_key == 1L) exp(d[["log_b[1]"]][i]) else NA
    for (r in seq_len(sd$n_zc)) {
      eta <- if (sd$K_det) sum(sd$Z_det[r, ] * bd) else 0
      ls  <- d$log_sigma[i] + eta + (if (sd$use_size_covar) d$beta_size[i] * (seq_len(K) - sd$s_centre) else 0)
      esw <- vapply(exp(ls), function(sg) integrate(function(x)
        if (sd$det_key == 1L) -expm1(-(sg / pmax(x, 1e-9))^b) else exp(-x^2 / (2 * sg^2)),
        0, sd$w, rel.tol = 1e-8)$value, 0)
      ref <- if (sd$use_size_covar) sum(pk * esw) else esw[1]
      out <- c(out, d[[paste0("esw_zc[", r, "]")]][i] / ref - 1)
    }
  }
  max(abs(out))
}

ok <- TRUE
for (spn in c("humpback", "pwsd")) {
  sd <- format_stan_data_visual(sim, spn, HSGP_M = c(6, 4), use_gp = 0L, use_bathy = FALSE)
  fit <- mod$sample(data = sd, chains = 2, parallel_chains = 2, iter_warmup = 300,
                    iter_sampling = 300, seed = 5, init = visual_init(sd), refresh = 0,
                    show_messages = FALSE, adapt_delta = 0.9)
  p <- sim$truth$sp_params[[spn]]
  truth <- c(log_sigma = log(p$sigma_det), beta_size = p$beta_size_truth,
             if (sd$det_key == 1L) c(b_shape = p$hr_b_truth),
             if (sd$K_det > 0) setNames(p$beta_det_truth, paste0("beta_det[", seq_len(sd$K_det), "]")))
  draws <- as_draws_matrix(fit$draws(names(truth)))
  lo <- apply(draws, 2, quantile, 0.005); hi <- apply(draws, 2, quantile, 0.995)
  inside <- truth >= lo & truth <= hi
  rel <- esw_rel_diff(fit, sd)
  cat(sprintf("\n%s (n = %d, n_zc = %d, divergences = %d)\n", spn, sd$n, sd$n_zc,
              sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)))
  print(data.frame(truth = truth, mean = colMeans(draws), lo99 = lo, hi99 = hi, inside = inside))
  cat(sprintf("max relative ESW difference vs integrate(): %.2g\n", rel))
  ok <- ok && all(inside) && rel < 1e-3
}
cat(if (ok) "\nPASS\n" else "\nFAIL\n")
quit(status = if (ok) 0L else 1L)
