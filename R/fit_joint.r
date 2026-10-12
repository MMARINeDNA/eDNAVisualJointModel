# =============================================================================
# fit_joint.r
#
# Fit stan/hsgp_joint.stan to a simulate_joint() dataset in one of three modes
# and summarise convergence + recovery against the simulated truth:
#
#   mode = "joint"   eDNA + visual likelihoods
#   mode = "edna"    eDNA only     (use_visual = 0)
#   mode = "visual"  visual only   (use_edna = 0)
#
# All three use the same model, priors, basis and initial values, so they
# differ only in the data that enter.
#
#   source("R/fit_joint.r")
#   sim <- simulate_joint(seed = 1)
#   res <- fit_joint(sim, "joint")
#   res$field; res$recovery
# =============================================================================

suppressMessages({
  library(cmdstanr)
  library(posterior)
})
source("R/functions_joint.R")
options(cmdstanr_warn_inits = FALSE)   # zero-size parameters are left out of inits

JOINT_STAN    <- "stan/hsgp_joint.stan"
JOINT_INCLUDE <- c("stan/include", file.path(HSGP4EDNA_DIR, "stan", "include"))

compile_joint <- function(...) cmdstan_model(JOINT_STAN, include_paths = JOINT_INCLUDE, ...)

joint_init <- function(sd) {
  S <- sd$S; L <- sd$L
  function() {
    init <- list(
      mu_sp = sd$prior_mu_sp_mu + c(rnorm(S, 0, 0.1), 2 + runif(1, 0, 0.1)),   # junk ~ 4
      gp_sigma = runif(S, 0.8, 1.2),
      gp_l_raw = matrix(runif(S * sd$D1, 0.2, 0.4), S, sd$D1),
      z_beta = matrix(0, S, sd$M), beta_bathy = matrix(0, S, sd$K_bathy),
      beta0_phi = 0, gamma0_phi = 2, gamma1_phi = 1,
      g_raw = matrix(0, sd$N * sd$use_bottle_noise, S),
      cv_bottle = rep(0.3, sd$use_bottle_noise * (if (sd$bottle_theta_shared == 1L) 1L else S)),
      g_junk_raw = rep(0, sd$N * sd$use_junk_noise), cv_junk = rep(0.3, sd$use_junk_noise),
      log_conv_free = sd$prior_log_conv_mu,
      # hazard-rate species start their scale at about w / 2 (as fit_visual.r)
      log_sigma = ifelse(sd$det_key == 1L, log(unlist(sd$w) / 2), sd$log_sigma_prior_mean) +
                  rnorm(L, 0, 0.1),
      beta_size = rnorm(L, 0, 0.005),
      log_b = (sd$log_b_prior_mean + rnorm(L, 0, 0.1))[sd$det_key == 1L],
      mu_s = sd$mu_s_prior_shape / sd$mu_s_prior_rate, phi_s = runif(L, 0.5, 2),
      mu_log_s = sd$mu_log_prior_mean + rnorm(L, 0, 0.1),
      sigma_log_s = sd$sigma_log_prior_shape / sd$sigma_log_prior_rate)
    # zero-size parameters must be left out of cmdstanr inits
    init[lengths(init) > 0]
  }
}

MODE_SWITCHES <- list(joint = c(1L, 1L), edna = c(1L, 0L), visual = c(0L, 1L))

fit_joint <- function(sim, mode = c("joint", "edna", "visual"), conv_free = character(0),
                      stan_data = NULL, mod = NULL,
                      chains = 3L, warmup = 800L, sample = 800L,
                      adapt_delta = 0.85, max_treedepth = 12L, seed = 42L,
                      refresh = 0L, output_dir = NULL, ...) {
  mode <- match.arg(mode)
  sw <- MODE_SWITCHES[[mode]]
  sd <- if (is.null(stan_data))
    format_stan_data_joint(sim, use_edna = sw[1], use_visual = sw[2], conv_free = conv_free, ...) else stan_data
  mod <- if (is.null(mod)) compile_joint() else mod
  set.seed(seed)   # initial values are drawn in R
  t0  <- Sys.time()
  fit <- mod$sample(data = sd, chains = chains, parallel_chains = chains,
                    iter_warmup = warmup, iter_sampling = sample,
                    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                    seed = seed, init = joint_init(sd), refresh = refresh,
                    show_messages = FALSE, output_dir = output_dir)
  runtime_min <- as.numeric(Sys.time() - t0, units = "mins")
  c(list(fit = fit, stan_data = sd, mode = mode, runtime_min = runtime_min),
    summarise_joint(fit, sim, sd, mode))
}

# True log density at every field location (stations, then segments) and the
# posterior at the same locations (each station via its first sample row).
summarise_joint <- function(fit, sim, sd, mode) {
  sp   <- attr(sd, "species"); S <- sd$S
  locs <- sim$field$locations
  truth <- sim$field$log_lambda_si                          # stations then segments
  smp  <- sim$edna$design$samples
  st_row <- match(seq_len(sum(locs$source == "station")), smp$station)  # first sample per station

  dm <- function(v) fit$draws(v, format = "draws_matrix")
  ll_st  <- dm("log_lambda");     ll_sg <- dm("log_lambda_seg")
  col <- function(m, nr, rows, s) m[, (s - 1) * nr + rows, drop = FALSE]   # column-major [row, s]
  post_mean <- matrix(NA_real_, nrow(truth), S, dimnames = list(NULL, sp))
  field <- do.call(rbind, lapply(seq_len(S), function(s) {
    draws <- cbind(col(ll_st, sd$N, st_row, s), col(ll_sg, sd$n_seg, seq_len(sd$n_seg), s))
    pm <- colMeans(draws)
    post_mean[, s] <<- pm
    q  <- apply(draws, 2, quantile, c(0.025, 0.975))
    tr <- truth[, s]
    is_st <- locs$source == "station"
    data.frame(mode = mode, species = sp[s],
               R2 = cor(pm, tr)^2, R2_stations = cor(pm[is_st], tr[is_st])^2,
               R2_segments = cor(pm[!is_st], tr[!is_st])^2,
               rmse = sqrt(mean((pm - tr)^2)),
               ci_width = mean(q[2, ] - q[1, ]),
               ci_cover = mean(tr >= q[1, ] & tr <= q[2, ]))
  }))

  gp <- sim$edna$truth$gp_params
  vars <- c(sprintf("gp_l[%d,1]", 1:S), sprintf("gp_l[%d,2]", 1:S),
            sprintf("gp_sigma[%d]", 1:S), sprintf("mu_sp[%d]", 1:S))
  truth_v <- c(vapply(gp, `[[`, 0, "lx"), vapply(gp, `[[`, 0, "ly"),
               vapply(gp, `[[`, 0, "sigma"), vapply(gp, `[[`, 0, "mu"))
  if (sd$n_conv_free > 0) {
    vars <- c(vars, sprintf("log_conv[%d]", sd$conv_free_idx))
    truth_v <- c(truth_v, sd$log_conv_factor[sd$conv_free_idx])
  }
  if (sd$use_visual == 1L) {
    lt <- sim$meta$lt_species; p <- sim$visual$truth$sp_params
    vars <- c(vars, sprintf("sigma_c[%d]", seq_len(sd$L)))
    truth_v <- c(truth_v, vapply(lt, function(x) p[[x]]$sigma_det, 0))
  }
  noise <- c(if (sd$use_bottle_noise == 1L) "theta_bottle", if (sd$use_junk_noise == 1L) "theta_junk")
  qs <- summarise_draws(fit$draws(c(vars, noise)), median = stats::median,
                        q025 = ~as.numeric(quantile(.x, 0.025)),
                        q975 = ~as.numeric(quantile(.x, 0.975)),
                        rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  recovery <- merge(data.frame(variable = vars, truth = truth_v), as.data.frame(qs),
                    by = "variable", all.y = TRUE, sort = FALSE)
  # Shared bottle-noise theta: simulate_edna_obs()'s theta_bottle (10). theta_junk
  # has no exact truth (the junk draw is not a gamma multiplier).
  recovery$truth[recovery$variable == "theta_bottle[1]"] <- 10
  recovery$covered <- with(recovery, ifelse(is.na(truth), NA, truth >= q025 & truth <= q975))
  recovery <- cbind(mode = mode, recovery)

  diag <- fit$diagnostic_summary(quiet = TRUE)
  hyper_rhat <- fit$summary(c("mu_sp", "gp_sigma", "gp_l_raw", "beta_bathy",
                              if (sd$use_visual == 1L) c("log_sigma", "mu_s", "mu_log_s"),
                              if (sd$n_conv_free > 0) "log_conv_free", noise), "rhat")$rhat
  list(field = field, recovery = recovery, post_mean = post_mean,
       divergences = sum(diag$num_divergent), treedepth_hits = sum(diag$num_max_treedepth),
       ebfmi = diag$ebfmi, max_rhat = max(hyper_rhat, na.rm = TRUE))
}
