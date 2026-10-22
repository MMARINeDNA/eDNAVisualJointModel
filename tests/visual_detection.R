# =============================================================================
# tests/visual_detection.R
#
# Checks the detection-function options of stan/hsgp_visual.stan on simulated
# data. The spatial field is ON (small basis): simulate_visual(sea_state = TRUE)
# makes Beaufort spatially correlated with density, so a non-spatial fit lets
# the Beaufort effect absorb density variation (confounding, not a model bug).
#
#   humpback-like   hazard-rate key + group-size covariate (Beaufort and
#                   visibility effects 0)
#   pwsd-like       half-normal key + group size + Beaufort + visibility
#                   (sea_state = TRUE, visibility = TRUE)
#
# For each: (1) the true detection parameters lie inside the 99% posterior
# interval, and (2) Stan's per-cell (Beaufort x visibility) ESW agrees with R's integrate()
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
                                            use_size_covar = 1L, beta_size_truth = 0.1,
                                            beta_bft_truth = 0, beta_vis_truth = 0))
sp$pwsd <- modifyList(sp$pwsd, list(det_key = 0L, sigma_det = 1.5, beta_size_truth = 0.005,
                                    beta_bft_truth = -0.15, beta_vis_truth = 0.2))
gp <- default_gp_params(); gp$pwsd$mu <- log(0.3); gp$humpback$mu <- log(0.02)
sim <- simulate_visual(seed = 11, n_transects = 25, seg_length = 10, w = 6.5,
                       gp_params = gp, sp_params = sp, sea_state = TRUE, visibility = TRUE)

esw_rel_diff <- function(fit, sd) {
  vars <- c("log_sigma", "beta_size", if (sd$use_bft_covar == 1L) "beta_bft",
            if (sd$use_vis_covar == 1L) "beta_vis",
            if (sd$det_key == 1L) "log_b", "pk", "esw_bft")
  d <- as.data.frame(fit$draws(vars, format = "draws_df"))
  out <- numeric(0)
  for (i in 1:5) {
    pk <- as.numeric(d[i, grep("^pk\\[", names(d))]); K <- length(pk)
    b  <- if (sd$det_key == 1L) exp(d[["log_b[1]"]][i]) else NA
    bb <- if (sd$use_bft_covar == 1L) d[["beta_bft[1]"]][i] else 0
    bv <- if (sd$use_vis_covar == 1L) d[["beta_vis[1]"]][i] else 0
    for (lv in seq_len(sd$n_bft)) {
      ls  <- d$log_sigma[i] + bb * sd$bft_level[lv] + bv * sd$vis_level[lv] +
        (if (sd$use_size_covar) d$beta_size[i] * (seq_len(K) - sd$s_centre) else 0)
      esw <- vapply(exp(ls), function(sg) integrate(function(x)
        if (sd$det_key == 1L) -expm1(-(sg / pmax(x, 1e-9))^b) else exp(-x^2 / (2 * sg^2)),
        0, sd$w, rel.tol = 1e-8)$value, 0)
      ref <- if (sd$use_size_covar) sum(pk * esw) else esw[1]
      out <- c(out, d[[paste0("esw_bft[", lv, "]")]][i] / ref - 1)
    }
  }
  max(abs(out))
}

ok <- TRUE
for (spn in c("humpback", "pwsd")) {
  sd <- format_stan_data_visual(sim, spn, HSGP_M = c(6, 4), use_gp = 1L, use_bathy = FALSE,
                                use_bft = TRUE, use_vis = spn == "pwsd")
  fit <- mod$sample(data = sd, chains = 2, parallel_chains = 2, iter_warmup = 300,
                    iter_sampling = 300, seed = 5, init = visual_init(sd), refresh = 0,
                    show_messages = FALSE, adapt_delta = 0.9)
  p <- sim$truth$sp_params[[spn]]
  truth <- c(log_sigma = log(p$sigma_det), beta_size = p$beta_size_truth,
             `beta_bft[1]` = p$beta_bft_truth,
             if (sd$use_vis_covar == 1L) c(`beta_vis[1]` = p$beta_vis_truth),
             if (sd$det_key == 1L) c(b_shape = p$hr_b_truth))
  draws <- as_draws_matrix(fit$draws(names(truth)))
  lo <- apply(draws, 2, quantile, 0.005); hi <- apply(draws, 2, quantile, 0.995)
  inside <- truth >= lo & truth <= hi
  rel <- esw_rel_diff(fit, sd)
  cat(sprintf("\n%s (n = %d, n_bft = %d, divergences = %d)\n", spn, sd$n, sd$n_bft,
              sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)))
  print(data.frame(truth = truth, mean = colMeans(draws), lo99 = lo, hi99 = hi, inside = inside))
  cat(sprintf("max relative ESW difference vs integrate(): %.2g\n", rel))
  ok <- ok && all(inside) && rel < 1e-3
}
cat(if (ok) "\nPASS\n" else "\nFAIL\n")
quit(status = if (ok) 0L else 1L)
