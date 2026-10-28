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

compile_joint <- function(...) stan_model_fresh(JOINT_STAN, JOINT_INCLUDE, ...)

# -----------------------------------------------------------------------------
# Data-based starting field for the eDNA stations.
#
# From a flat starting field (z = 0) the sampler can get stuck where the data
# imply extreme densities: at hot spots it settles on near-zero copies with
# inflated noise (joint simulator seed 1: hake R^2 0.19, theta_bottle 1.6
# against 10, log density ~2000 below the true region). Starting the field
# from the data fixes it (hake R^2 0.978, theta_bottle 8.9; the same as
# starting at the truth). ROADMAP_2026-10-06, 6.0.
#
# Per station and species, a target log density from the data only:
#   - metabarcoding: log(reads_s / (reads_junk + 0.5)) + the junk's starting
#     log copies, averaged over the replicates with reads of species s;
#   - qPCR (species 1): log aliquot copies from Ct through the standard curve
#     (with estimated curves, each well's plate curve fitted to the standards
#     alone), minus the well's dilution, averaged over detected replicates
#     (replaces the MB value);
#   - minus the conversion factor, the bottle's log volume and the
#     water-column offset (conversion factors being estimated use their
#     starting value).
# As HSGP4eDNA's field_start_from_data(); with no dilutions and one volume
# it is the start before v1.8 exactly.
# Stations without a value get the species' minimum minus 1. The targets are smoothed
# onto the HSGP basis by the GP posterior mean of z (ridge, noise SD 0.5) at
# the prior's median length-scale, giving the starting mu, gp_sigma, gp_l_raw
# and z_beta. Uses the Stan data only, so it applies to real data too.
# -----------------------------------------------------------------------------
edna_field_start <- function(sd, mu_junk, log_conv_free = sd$prior_log_conv_mu) {
  S <- sd$S; N <- sd$N; M <- sd$M; D <- sd$D1; L <- sd$L_hsgp; IDX <- sd$INDICES
  PHI <- sapply(seq_len(M), function(m) Reduce(`*`, lapply(seq_len(D), function(i)
    1 / sqrt(L[i]) * sin(IDX[m, i] * pi * (sd$coords[, i] + L[i]) / (2 * L[i])))))
  sqrt_spd <- function(alpha, rho) apply(IDX, 1, function(mm)
    sqrt(alpha^2 * sqrt(2 * pi)^D * prod(rho) * exp(-0.5 * sum(rho^2 * (mm * pi / (2 * L))^2))))
  rho0 <- sd$prior_gp_l_scale / qgamma(0.5, sd$prior_gp_l_shape)   # inverse-gamma median

  log_conv <- sd$log_conv_factor
  if (sd$n_conv_free > 0) log_conv[sd$conv_free_idx] <- log_conv_free
  lvf <- log(sd$vol_aliquot / 100)
  log_vol <- sd$log_vol_bottle
  dens <- function(log_bottle, s, i)
    log_bottle - log_conv[s] - log_vol[i] - sd$log_zsample_effect[i, s]
  target <- matrix(NA_real_, N, S)
  log_junk <- mu_junk + log_conv[S + 1] + log_vol
  for (s in 1:S) {
    pos <- sd$mb_reads[, s] > 0
    if (!any(pos)) next
    v <- tapply(log(sd$mb_reads[pos, s]) - log(sd$mb_reads[pos, S + 1] + 0.5), sd$mb_sample_idx[pos], mean)
    i <- as.integer(names(v)); target[i, s] <- dens(log_junk[i] + v, s, i)
  }
  det <- sd$qpcr_detect == 1
  if (any(det)) {
    if (sd$qpcr_curve_model == 1L) {
      cs <- qpcr_curve_start(sd); pl <- sd$qpcr_plate_idx[det]
      log_well <- (cs$alpha_ct[pl] - sd$qpcr_ct[det]) / cs$beta_ct[pl]
    } else log_well <- (sd$alpha_ct - sd$qpcr_ct[det]) / sd$beta_ct
    v <- tapply(log_well - lvf - sd$qpcr_log_dilution[det], sd$qpcr_sample_idx[det], mean)
    i <- as.integer(names(v)); target[i, 1] <- dens(v, 1, i)
  }
  out <- list(mu = numeric(S), sigma = numeric(S), rho = rho0, z = matrix(0, S, M))
  for (s in 1:S) {
    y <- target[, s]
    y[is.na(y)] <- if (all(is.na(y))) sd$prior_mu_sp_mu[s] else min(y, na.rm = TRUE) - 1
    out$mu[s] <- mean(y); out$sigma[s] <- max(sd(y), 0.5)
    A <- sweep(PHI, 2, sqrt_spd(out$sigma[s], rho0), `*`)
    out$z[s, ] <- t(A) %*% solve(tcrossprod(A) + diag(0.25, N), y - out$mu[s])
  }
  out
}

joint_init <- function(sd) {
  S <- sd$S; L <- sd$L
  mu_junk <- sd$prior_mu_sp_mu[S + 1] + 2.05                               # junk ~ 4
  # data-based starting field when the eDNA data are in the likelihood (a
  # visual-only fit must not start from the eDNA data)
  fs <- if (sd$use_edna == 1L) edna_field_start(sd, mu_junk) else NULL
  cs <- if (sd$qpcr_curve_model == 1L) qpcr_curve_start(sd) else NULL
  function() {
    init <- list(
      mu_sp = c(if (is.null(fs)) sd$prior_mu_sp_mu[1:S] else fs$mu, mu_junk) +
              c(rnorm(S, 0, 0.1), runif(1, -0.05, 0.05)),
      gp_sigma = (if (is.null(fs)) 1 else fs$sigma) * runif(S, 0.8, 1.2),
      gp_l_raw = if (is.null(fs)) matrix(runif(S * sd$D1, 0.2, 0.4), S, sd$D1)
                 else matrix(rep(fs$rho, each = S), S, sd$D1) * runif(S * sd$D1, 0.9, 1.1),
      z_beta = if (is.null(fs)) matrix(0, S, sd$M) else fs$z + rnorm(S * sd$M, 0, 0.05),
      beta_bathy = matrix(0, S, sd$K_bathy),
      beta0_phi = 0, gamma0_phi = 2, gamma1_phi = 1,
      g_raw = matrix(0, sd$N * sd$use_bottle_noise, S),
      cv_bottle = rep(0.3, sd$use_bottle_noise * (if (sd$bottle_theta_shared == 1L) 1L else S)),
      g_junk_raw = rep(0, sd$N * sd$use_junk_noise), cv_junk = rep(0.3, sd$use_junk_noise),
      log_conv_free = sd$prior_log_conv_mu,
      beta_zdepth = matrix(0, S, sd$K_zdepth),
      # hazard-rate species start their scale at about w / 2 (as fit_visual.r)
      log_sigma = ifelse(sd$det_key == 1L, log(unlist(sd$w) / 2), sd$log_sigma_prior_mean) +
                  rnorm(L, 0, 0.1),
      beta_size = rnorm(L, 0, 0.005),
      log_b = pmax(sd$log_b_prior_mean + rnorm(L, 0, 0.1), log(sd$b_min) + 0.05)[sd$det_key == 1L],
      # sea-state slopes start near 0 (draws nothing when no species uses it)
      beta_bft = array(rnorm(sum(sd$use_bft_covar), 0, 0.02)),
      beta_vis = array(rnorm(sum(sd$use_vis_covar), 0, 0.02)),
      inv_sqrt_phi_nb = array(runif(sum(sd$use_nb), 0.3, 0.7)),
      mu_s = sd$mu_s_prior_shape / sd$mu_s_prior_rate, phi_s = runif(L, 0.5, 2),
      mu_log_s = sd$mu_log_prior_mean + rnorm(L, 0, 0.1),
      sigma_log_s = sd$sigma_log_prior_shape / sd$sigma_log_prior_rate)
    # per-plate qPCR curves start from the standards alone (as HSGP4eDNA's edna_init())
    if (!is.null(cs)) init <- c(init, list(
      alpha_ct_plate = cs$alpha_ct + rnorm(sd$N_plate, 0, 0.05),
      beta_ct_plate = pmax(cs$beta_ct + rnorm(sd$N_plate, 0, 0.002), sd$beta_ct_min + 1e-4),
      gamma0_ct_est = array(0, 1), gamma1_ct_est = array(-0.2, 1),
      sigma0_ct_est = array(0.3, 1), logit_kappa_est = array(2, 1)))
    # zero-size parameters must be left out of cmdstanr inits
    init[lengths(init) > 0]
  }
}

# "prior": neither likelihood (the prior predictive; the eDNA noise terms are
# kept, drawn from their priors)
MODE_SWITCHES <- list(joint = c(1L, 1L), edna = c(1L, 0L), visual = c(0L, 1L), prior = c(0L, 0L))

fit_joint <- function(sim, mode = c("joint", "edna", "visual", "prior"), conv_free = character(0),
                      stan_data = NULL, mod = NULL,
                      chains = 4L, warmup = 800L, sample = 800L,
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
                    show_messages = refresh > 0, output_dir = output_dir)
  runtime_min <- as.numeric(Sys.time() - t0, units = "mins")
  c(list(fit = fit, stan_data = sd, mode = mode, runtime_min = runtime_min),
    if (is.null(sim$field)) summarise_joint_real(fit, sim, sd, mode) else summarise_joint(fit, sim, sd, mode))
}

# Hyperparameters whose R-hat the summaries report (the eDNA options' when on)
joint_hyper_vars <- function(sd) {
  noise <- c(if (sd$use_bottle_noise == 1L) "theta_bottle", if (sd$use_junk_noise == 1L) "theta_junk")
  c("mu_sp", "gp_sigma", "gp_l_raw", "beta_bathy",
    if (sd$use_visual == 1L) c("log_sigma", "mu_s", "mu_log_s"),
    if (sd$use_visual == 1L && sum(sd$use_bft_covar) > 0) "beta_bft",
    if (sd$use_visual == 1L && sum(sd$use_nb) > 0) "inv_sqrt_phi_nb",
    if (sd$use_visual == 1L && sum(sd$use_vis_covar) > 0) "beta_vis",
    if (sd$n_conv_free > 0) "log_conv_free",
    if (sd$use_edna == 1L && sd$qpcr_curve_model == 1L)
      c("alpha_ct_plate", "beta_ct_plate", "gamma0_ct_est", "gamma1_ct_est", "sigma0_ct_est", "logit_kappa_est"),
    if (sd$use_edna == 1L && sd$K_zdepth > 0) "beta_zdepth", noise)
}

# Real data (no truth): posterior mean and 90% interval of the field at every
# station and segment, the main parameters, the length-scale floor check and
# diagnostics.
summarise_joint_real <- function(fit, sim, sd, mode) {
  sp <- attr(sd, "species"); S <- sd$S
  smp <- sim$edna$design$samples
  st_row <- match(sort(unique(smp$station)), smp$station)        # first sample per station
  dm <- function(v) fit$draws(v, format = "draws_matrix")
  col <- function(m, nr, rows, s) m[, (s - 1) * nr + rows, drop = FALSE]
  ll_st <- dm("log_lambda"); ll_sg <- dm("log_lambda_seg")
  field <- lapply(setNames(seq_len(S), sp), function(s) {
    st <- col(ll_st, sd$N, st_row, s); sg <- col(ll_sg, sd$n_seg, seq_len(sd$n_seg), s)
    q <- function(m) t(apply(m, 2, quantile, c(0.05, 0.95)))
    list(station = data.frame(station = smp$station[st_row], X = smp$X[st_row], Y = smp$Y[st_row],
                              mean = colMeans(st), q = q(st)),
         segment = data.frame(seg_id = seq_len(sd$n_seg), mean = colMeans(sg), q = q(sg)))
  })
  vars <- c("mu_sp", "gp_sigma", "gp_l", "beta_bathy",
            if (sd$n_conv_free > 0) "log_conv",
            if (sd$use_visual == 1L) c("sigma_c", "b_shape", "p_det", "mean_group_size"),
            if (sd$use_visual == 1L && sum(sd$use_bft_covar) > 0) "beta_bft",
            if (sd$use_visual == 1L && sum(sd$use_vis_covar) > 0) "beta_vis",
            if (sd$use_visual == 1L && sum(sd$use_nb) > 0) "phi_nb",
            if (sd$use_bottle_noise == 1L) "theta_bottle", if (sd$use_junk_noise == 1L) "theta_junk",
            if (sd$use_edna == 1L && sd$qpcr_curve_model == 1L)
              c("alpha_ct_plate", "beta_ct_plate", "kappa_qpcr", "gamma0_ct_est", "gamma1_ct_est", "sigma0_ct_est"),
            if (sd$use_edna == 1L && sd$K_zdepth > 0) "beta_zdepth",
            if (sd$use_edna == 1L) c("beta0_phi", "gamma0_phi", "gamma1_phi"))
  params <- cbind(mode = mode, as.data.frame(summarise_draws(fit$draws(vars), "mean", median = stats::median,
                  q05 = ~as.numeric(quantile(.x, 0.05)), q95 = ~as.numeric(quantile(.x, 0.95)),
                  rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)))
  # each chain's mean of every parameter and lp__ (variables x chains), to see
  # which chain sits apart when R-hat is high
  chain_means <- apply(posterior::as_draws_array(fit$draws(c("lp__", vars))), c(3, 2), mean)
  gp_prior <- attr(sd, "gp_prior")
  floor <- do.call(rbind, lapply(seq_len(S), function(s) cbind(mode = mode, species = sp[s],
             floor_check(fit$draws(sprintf("gp_l[%d,%d]", s, 1:2), format = "draws_matrix"), gp_prior))))
  diag <- fit$diagnostic_summary(quiet = TRUE)
  hyper_rhat <- fit$summary(joint_hyper_vars(sd), "rhat")$rhat
  list(field = field, params = params, floor = floor, chain_means = chain_means,
       grid = summarise_grid(fit, sd), ppc = ppc_joint(fit, sim, sd),
       divergences = sum(diag$num_divergent), treedepth_hits = sum(diag$num_max_treedepth),
       ebfmi = diag$ebfmi, max_rhat = max(hyper_rhat, na.rm = TRUE))
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
    # Beaufort slopes (beta_bft is indexed over the species that use it); the
    # truth is the simulator's when the segments carry sea state, else 0
    bft_lt <- lt[sd$use_bft_covar == 1L]
    if (length(bft_lt)) {
      vars <- c(vars, sprintf("beta_bft[%d]", seq_along(bft_lt)))
      truth_v <- c(truth_v, vapply(bft_lt, function(x)
        if (isTRUE(sim$visual$meta$sea_state)) p[[x]]$beta_bft_truth else 0, 0))
    }
    # Visibility slopes (indexed over the species that use it)
    vis_lt <- lt[sd$use_vis_covar == 1L]
    if (length(vis_lt)) {
      vars <- c(vars, sprintf("beta_vis[%d]", seq_along(vis_lt)))
      truth_v <- c(truth_v, vapply(vis_lt, function(x)
        if (isTRUE(sim$visual$meta$visibility)) p[[x]]$beta_vis_truth else 0, 0))
    }
    # NB sizes (indexed over the NB species); truth NA for Poisson sightings
    nb_lt <- lt[sd$use_nb == 1L]
    if (length(nb_lt)) {
      vars <- c(vars, sprintf("phi_nb[%d]", seq_along(nb_lt)))
      truth_v <- c(truth_v, vapply(nb_lt, function(x) {
        k <- p[[x]]$nb_size_truth; if (is.null(k) || !is.finite(k)) NA_real_ else k }, 0))
    }
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
  hyper_rhat <- fit$summary(joint_hyper_vars(sd), "rhat")$rhat
  list(field = field, recovery = recovery, post_mean = post_mean,
       grid = summarise_grid(fit, sd), ppc = ppc_joint(fit, sim, sd),
       divergences = sum(diag$num_divergent), treedepth_hits = sum(diag$num_max_treedepth),
       ebfmi = diag$ebfmi, max_rhat = max(hyper_rhat, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# Field on the prediction grid (log_lambda_pred): posterior mean, SD and 90%
# interval per point, one data frame per species. NULL without a grid.
# -----------------------------------------------------------------------------
summarise_grid <- function(fit, sd) {
  if (sd$N_pred == 0L) return(NULL)
  sp <- attr(sd, "species")
  setNames(lapply(seq_len(sd$S), function(s) {
    m <- fit$draws(sprintf("log_lambda_pred[%d,%d]", seq_len(sd$N_pred), s), format = "draws_matrix")
    q <- apply(m, 2, quantile, c(0.05, 0.95))
    data.frame(point = seq_len(sd$N_pred), mean = colMeans(m), sd = apply(m, 2, stats::sd),
               q05 = q[1, ], q95 = q[2, ])
  }), sp)
}

# -----------------------------------------------------------------------------
# Predictive checks, computed from n_draws posterior (or prior) draws and kept
# as small summaries, since the fits themselves are not saved.
#
# Visual (visual-only, joint and prior fits), per line-transect species:
#   counts   segment counts replicated from the fitted means (Poisson or NB):
#            total, share of zero segments, variance / mean across segments,
#            maximum; per draw (rep) and observed (obs)
#   bands    detected groups per 100 km band of Y (along the coast): observed
#            against the replicated median and 90% interval
#   dist     the predictive density of detected perpendicular distance
#            (det_pdf) against the observed distances
#   size     detected group sizes by size class: predictive (det_size)
#            against observed shares
# eDNA (eDNA-only, joint and prior fits):
#   mb       share of metabarcoding replicates with reads, per species and
#            sampling depth: P(reads > 0) = (1 - exp(-aliquot copies)) x
#            (1 - BetaBinomial(0 | total)), as the likelihood
#   qpcr     hake qPCR detection rate per sampling depth
#   ct       per detected well: observed Ct, predictive mean and 90% interval
#            (the Ct noise as the model's, with ct_sd_cond)
# -----------------------------------------------------------------------------
ppc_joint <- function(fit, sim, sd, n_draws = 200L) {
  nd_all <- posterior::ndraws(fit$draws("lp__"))
  sel <- unique(round(seq(1, nd_all, length.out = min(n_draws, nd_all))))
  prior <- sd$use_edna == 0L && sd$use_visual == 0L
  list(visual = if (sd$L > 0 && (sd$use_visual == 1L || prior)) ppc_visual(fit, sim, sd, sel),
       edna   = if (sd$use_edna == 1L || prior) ppc_edna(fit, sim, sd, sel))
}

# group-size bins of the model (sizes 1..30, then 10% log bins) and classes
ksz_bin_r <- function(k) ifelse(k <= 30, k, 31 + floor(log(k) / log(1.1)) - floor(log(31) / log(1.1)))
size_classes <- function(S_max, edges = c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000)) {
  k <- seq_len(S_max); b <- ksz_bin_r(k)
  first <- tapply(k, b, min)                       # first size in each bin
  edges <- unique(first[findInterval(edges[edges <= S_max], first)])   # moved to bin starts
  cls <- findInterval(first, edges)
  last <- c(edges[-1] - 1, S_max)
  list(bin_class = cls, n_bin = length(first), edges = edges,
       label = ifelse(edges == last, as.character(edges), sprintf("%d-%d", edges, last)),
       size_class = function(s) findInterval(first[ksz_bin_r(s)], edges))
}

ppc_visual <- function(fit, sim, sd, sel) {
  dm <- function(v) fit$draws(v, format = "draws_matrix")[sel, , drop = FALSE]
  seg <- sim$visual$design$segments
  band <- floor(seg$Y / 100)
  mgs <- dm("mean_group_size"); nd <- length(sel)
  lt <- sim$meta$lt_species
  det_end <- cumsum(sd$n_det)
  setNames(lapply(seq_len(sd$L), function(l) {
    s <- sd$lt_species_idx[l]
    ll  <- dm(sprintf("log_lambda_seg[%d,%d]", seq_len(sd$n_seg), s))
    esw <- dm(sprintf("esw_bft[%d,%d]", l, seq_len(sd$n_bft)))[, sd$seg_bft, drop = FALSE]
    mu  <- exp(unclass(ll) - log(as.numeric(mgs[, l])) + log(unclass(esw)) + rep(log(2 * as.numeric(sd$seg_l)), each = nd))
    yrep <- if (sd$use_nb[l] == 1L) {
      phi <- as.numeric(dm(sprintf("phi_nb[%d]", cumsum(sd$use_nb)[l])))
      matrix(rnbinom(length(mu), size = rep(phi, ncol(mu)), mu = mu), nd)
    } else matrix(rpois(length(mu), mu), nd)
    y <- sd$seg_count[l, ]
    stat <- function(m) data.frame(total = rowSums(m), zero = rowMeans(m == 0),
                                   disp = apply(m, 1, var) / pmax(rowMeans(m), 1e-12), max = apply(m, 1, max))
    br <- t(rowsum(t(yrep), band)); bo <- tapply(y, band, sum)
    bq <- apply(br, 2, quantile, c(0.05, 0.5, 0.95))
    out <- list(counts = list(rep = stat(yrep), obs = stat(matrix(y, 1))),
                bands = data.frame(Y_km = (as.numeric(colnames(br)) + 0.5) * 100, obs = as.numeric(bo),
                                   lo = bq[1, ], med = bq[2, ], hi = bq[3, ]))
    if (sd$use_ppc == 1L) {
      i <- (det_end[l] - sd$n_det[l] + 1):det_end[l]
      pdf <- dm(sprintf("det_pdf[%d,%d]", l, 1:40)); pq <- apply(pdf, 2, quantile, c(0.05, 0.5, 0.95))
      out$dist <- list(pred = data.frame(x = (1:40 - 0.5) * sd$w[l] / 40, lo = pq[1, ], med = pq[2, ], hi = pq[3, ]),
                       obs = sd$det_x[i], w = sd$w[l])
      sc <- size_classes(sd$S_max[l])
      pb <- dm(sprintf("det_size[%d,%d]", l, seq_len(sc$n_bin)))
      pc <- t(rowsum(t(pb), sc$bin_class)); cq <- apply(pc, 2, quantile, c(0.05, 0.5, 0.95))
      oc <- tabulate(sc$size_class(sd$det_s[i]), length(sc$edges)) / length(i)
      out$size <- data.frame(class = factor(sc$label, sc$label), obs = oc,
                             lo = cq[1, ], med = cq[2, ], hi = cq[3, ], n_obs = length(i))
    }
    out
  }), lt)
}

ppc_edna <- function(fit, sim, sd, sel) {
  dm <- function(v) fit$draws(v, format = "draws_matrix")[sel, , drop = FALSE]
  nd <- length(sel); S <- sd$S; K <- S + 1; N <- sd$N
  le <- array(dm(sprintf("log_lambda_edna[%d,%d]", rep(seq_len(N), K), rep(seq_len(K), each = N))), c(nd, N, K))
  smp <- sim$edna$design$samples
  depth <- if (!is.null(smp$Z_sample)) smp$Z_sample else rep(0, N)
  sp <- attr(sd, "species")
  lvf <- log(sd$vol_aliquot / 100)
  # log of the expected copies given presence, lam / (1 - exp(-lam)), stable for small lam
  log_k <- function(x) { lam <- exp(x); ifelse(lam < 1e-8, lam / 2, x - log(-expm1(-lam))) }
  lse <- function(m) { mx <- apply(m, 1, max); mx + log(rowSums(exp(m - mx))) }
  b0 <- as.numeric(dm("beta0_phi")); g0 <- as.numeric(dm("gamma0_phi")); g1 <- as.numeric(dm("gamma1_phi"))
  ri <- sd$mb_sample_idx; R <- length(ri)
  mb_off <- lvf + as.numeric(sd$mb_log_dilution); mb_tot <- as.numeric(sd$mb_total)
  ppos <- array(NA_real_, c(nd, R, S))
  for (d in seq_len(nd)) {
    LL <- matrix(le[d, ri, ], R, K) + mb_off                                      # R x K
    lK <- log_k(LL); lnu <- lK - lK[, K]; lpi <- lnu - lse(lnu)
    lphi <- b0[d] + lse(LL) + pmax(g0[d] - g1[d] * lK, 0)
    a <- exp(lpi + lphi); b <- exp(log(-expm1(pmin(lpi, -1e-300))) + lphi)
    p_bb0 <- exp(lbeta(a, mb_tot + b) - lbeta(a, b))
    ppos[d, , ] <- (-expm1(-exp(LL[, 1:S]))) * (1 - p_bb0[, 1:S])
  }
  dz <- depth[ri]
  mb <- do.call(rbind, lapply(seq_len(S), function(s) do.call(rbind, lapply(sort(unique(dz)), function(z) {
    r <- dz == z; pr <- rowMeans(ppos[, r, s, drop = FALSE])
    data.frame(species = sp[s], depth = z, n = sum(r), obs = mean(sd$mb_reads[r, s] > 0),
               lo = unname(quantile(pr, 0.05)), med = median(pr), hi = unname(quantile(pr, 0.95)))
  }))))
  out <- list(mb = mb)
  if (sd$N_qpcr_long > 0) {
    qi <- sd$qpcr_sample_idx
    lw <- matrix(le[, qi, 1], nd) + rep(lvf + as.numeric(sd$qpcr_log_dilution), each = nd)  # nd x wells: log copies in the well
    curve <- sd$qpcr_curve_model == 1L
    kap <- if (curve) as.numeric(dm("kappa_qpcr[1]")) else rep(sd$kappa, nd)
    pdet <- -expm1(-kap * exp(lw))
    dq <- depth[qi]
    out$qpcr <- do.call(rbind, lapply(sort(unique(dq)), function(z) {
      r <- dq == z; pr <- rowMeans(pdet[, r, drop = FALSE])
      data.frame(depth = z, n = sum(r), obs = mean(sd$qpcr_detect[r]),
                 lo = unname(quantile(pr, 0.05)), med = median(pr), hi = unname(quantile(pr, 0.95)))
    }))
    pl <- sd$qpcr_plate_idx
    if (curve) {
      # one column per plate, then per well (draws subsetting drops repeated names)
      up <- sort(unique(pl)); w <- match(pl, up)
      al <- unclass(dm(sprintf("alpha_ct_plate[%d]", up)))[, w, drop = FALSE]
      be <- unclass(dm(sprintf("beta_ct_plate[%d]", up)))[, w, drop = FALSE]
      g0c <- as.numeric(dm("gamma0_ct_est[1]")); g1c <- as.numeric(dm("gamma1_ct_est[1]"))
      s0c <- as.numeric(dm("sigma0_ct_est[1]"))
    } else {
      al <- matrix(sd$alpha_ct, nd, length(qi)); be <- matrix(sd$beta_ct, nd, length(qi))
      g0c <- rep(sd$gamma0_ct, nd); g1c <- rep(sd$gamma1_ct, nd); s0c <- rep(sd$sigma0_ct, nd)
    }
    mct <- al - be * lw
    lsd <- if (sd$ct_sd_cond == 1L) log_k(lw) else lw
    sct <- sqrt(s0c^2 + exp(2 * (g0c + g1c * lsd)))
    ct_rep <- mct + sct * matrix(rnorm(length(mct)), nd)
    det <- sd$qpcr_detect == 1L
    cq <- apply(ct_rep[, det, drop = FALSE], 2, quantile, c(0.05, 0.95))
    out$ct <- data.frame(well = which(det), depth = dq[det], obs = sd$qpcr_ct[det],
                         pred = colMeans(mct[, det, drop = FALSE]), lo = cq[1, ], hi = cq[2, ])
  }
  out
}
