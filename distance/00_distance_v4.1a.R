# =============================================================================
# 00_distance_v4.1a.R
#
# Spatial extension of distance/00_distance_v4.1.R. The non-spatial v4.1
# pipeline fit a single scalar log_lambda_s per species; this pipeline
# adds a 3-D anisotropic HSGP latent field over segment locations,
# matching the spatial structure of stan/whale_edna_hsgp_v3.2.stan but
# with the line-transect observation model.
#
# Per species (humpback, PWSD), this script:
#
#   1. Re-uses the v4.1 simulation machinery: same domain, same E-W
#      systematic transects, same per-species GP / detection /
#      group-size truth.
#   2. Builds the HSGP data block (normalised coords, m_hsgp = (14, 8, 32),
#      L_hsgp = 1.5, gp_l priors centred on the per-species named scales).
#   3. Fits distance/distance_v4.1a.stan -- one fit per species.
#   4. Adds **spatial posterior-predictive checks** on top of the v4.1
#      diagnostic suite:
#        - per-segment posterior lambda_groups vs simulated truth
#        - per-segment posterior f (spatial component) vs simulated truth
#        - posterior-predictive segment counts vs observed counts
#        - GP hyperparameter recovery (gp_sigma, lx, ly, lz)
#
# Run from the project root:
#   Rscript distance/00_distance_v4.1a.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(cmdstanr)
  library(posterior)
  library(patchwork)
  library(viridis)
})

set.seed(7)

# -----------------------------------------------------------------------------
# 0. Helpers (verbatim from 00_distance_v4.1.R)
# -----------------------------------------------------------------------------

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

aniso_cov <- function(coords, sigma, lx, ly, lz, jitter = 1e-6) {
  d2 <- outer(coords[, 1], coords[, 1], function(a, b) ((a - b) / lx)^2) +
        outer(coords[, 2], coords[, 2], function(a, b) ((a - b) / ly)^2) +
        outer(coords[, 3], coords[, 3], function(a, b) ((a - b) / lz)^2)
  sigma * sigma * exp(-0.5 * d2) + diag(jitter, nrow(d2))
}

rnbinom_zt <- function(n, mu, size) {
  out <- integer(n)
  for (i in seq_len(n)) {
    repeat {
      r <- rnbinom(1, mu = mu, size = size)
      if (r > 0L) { out[i] <- r; break }
    }
  }
  out
}

# -----------------------------------------------------------------------------
# 1. Empirical group-size pools (same as v4.1)
# -----------------------------------------------------------------------------
load("distance/humpback_grpsz.RData")
load("distance/pwsd_grpsz.RData")
gs_pool_humpback <- as.integer(si.humpback)
gs_pool_pwsd     <- as.integer(si.pwsd)

OUTPUT_DIR <- "outputs/distance_v4.1a"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. Domain (matches v4.1)
# -----------------------------------------------------------------------------
X_km_max     <- 500
Y_km_max     <- 1270
Z_bathy_max  <- 3500
rotation_deg <- 25
rot          <- rotation_deg * pi / 180

X_prime_fn <- function(X_km, Y_km) {
  raw <- cos(rot) * (X_km / X_km_max) + sin(rot) * (Y_km / Y_km_max)
  300 * raw / (cos(rot) + sin(rot))
}

bathy_mean_fn <- function(X_km, Y_km) {
  xp      <- X_prime_fn(X_km, Y_km)
  abyssal <- 2500
  shelf   <- 80
  pmax(abyssal + (shelf - abyssal) / (1 + exp(-0.06 * (xp - 180))), 10)
}

# HSGP coord normalisation (same scheme as v4.1's
# scripts/03_format_stan_data_v4.1.r). All segment coords land in
# [-1, 1] under coord_centre / coord_scale derived from domain extents.
coord_centre <- c(X_km_max / 2, Y_km_max / 2, Z_bathy_max / 2)
coord_scale  <- c(X_km_max / 2, Y_km_max / 2, Z_bathy_max / 2)

HSGP_M <- c(14L, 8L, 32L)    # M_total = 3584 (same as v4.1)
HSGP_C <- c(1.5, 1.5, 1.5)   # boundary extension factor

# -----------------------------------------------------------------------------
# 3. Per-species params (mirrors 00_distance_v4.1.R but adds GP-prior
#    settings for the spatial submodel).
# -----------------------------------------------------------------------------

humpback_mean_gs <- mean(gs_pool_humpback)
pwsd_mean_gs     <- mean(gs_pool_pwsd)

pwsd_log_mu0    <- log(pwsd_mean_gs) -
                   0.5 * log(1 + (sd(gs_pool_pwsd) / pwsd_mean_gs)^2)
pwsd_log_sigma0 <- sqrt(log(1 + (sd(gs_pool_pwsd) / pwsd_mean_gs)^2))

species_params <- list(

  humpback = list(
    common_name      = "Humpback whale",
    sigma            = 1.0,
    lx               = 50,
    ly               = 300,
    lz               = 100,
    mu               = log(0.004),
    mean_group_size  = humpback_mean_gs,
    sim_group_dist   = "empirical",
    group_size_pool  = gs_pool_humpback,

    # Detection function
    sigma_det        = 2.5,
    use_size_covar   = 0L,
    beta_size_truth  = 0.0,
    s_centre         = humpback_mean_gs,

    # Group-size sub-model: zero-truncated NB
    model_group_dist = 0L,

    # Per-species priors for the LT side (same as v4.1)
    log_sigma_prior_mean    = log(2.5),
    log_sigma_prior_sd      = 0.6,
    beta_size_prior_mean    = 0.0,
    beta_size_prior_sd      = 0.05,
    mu_s_prior_shape        = 4,
    mu_s_prior_rate         = 4 / humpback_mean_gs,
    phi_s_prior_shape       = 1,
    phi_s_prior_rate        = 0.1,
    mu_log_prior_mean       = log(humpback_mean_gs),
    mu_log_prior_sd         = 1.0,
    sigma_log_prior_shape   = 2,
    sigma_log_prior_rate    = 2,
    S_max                   = 50L,

    # Spatial-field priors. mu_sp is the log GROUP density intercept
    # (matching v4.1's log_lambda_s convention). gp_l priors centred on
    # the per-species named length-scales with ~30% sd. gp_sigma uses
    # the v3.2 / v4.1 gamma(8, 4) (mode 1.75, mean 2.0).
    prior_mu_sp_mu          = log(0.004 / humpback_mean_gs),
    prior_mu_sp_sig         = 1.5,
    prior_gp_sigma_shape    = 8.0,
    prior_gp_sigma_rate     = 4.0,
    prior_gp_lx_mu          =  50,  prior_gp_lx_sig = 15,
    prior_gp_ly_mu          = 300,  prior_gp_ly_sig = 100,
    prior_gp_lz_mu          = 100,  prior_gp_lz_sig =  35
  ),

  pwsd = list(
    common_name      = "Pacific white-sided dolphin",
    sigma            = 1.3,
    lx               = 40,
    ly               = 300,
    lz               = 300,
    mu               = log(0.032),
    mean_group_size  = pwsd_mean_gs,
    sim_group_dist   = "empirical",
    group_size_pool  = gs_pool_pwsd,

    sigma_det        = 1.5,
    use_size_covar   = 1L,
    beta_size_truth  = 0.01,
    s_centre         = pwsd_mean_gs,

    model_group_dist = 1L,

    log_sigma_prior_mean    = log(1.5),
    log_sigma_prior_sd      = 0.6,
    beta_size_prior_mean    = 0.0,
    beta_size_prior_sd      = 0.05,
    mu_s_prior_shape        = 4,
    mu_s_prior_rate         = 4 / pwsd_mean_gs,
    phi_s_prior_shape       = 1,
    phi_s_prior_rate        = 0.1,
    mu_log_prior_mean       = pwsd_log_mu0,
    mu_log_prior_sd         = 0.5,
    sigma_log_prior_shape   = 4,
    sigma_log_prior_rate    = 4 / pwsd_log_sigma0,
    S_max                   = 1000L,

    prior_mu_sp_mu          = log(0.032 / pwsd_mean_gs),
    prior_mu_sp_sig         = 1.5,
    prior_gp_sigma_shape    = 8.0,
    prior_gp_sigma_rate     = 4.0,
    prior_gp_lx_mu          =  40,  prior_gp_lx_sig =  15,
    prior_gp_ly_mu          = 300,  prior_gp_ly_sig = 100,
    prior_gp_lz_mu          = 300,  prior_gp_lz_sig = 100
  )
)

w <- 5.0   # truncation distance, km

# -----------------------------------------------------------------------------
# 4. Survey design: 25 systematic E-W transects, 10 km segments
# -----------------------------------------------------------------------------
n_transects        <- 25
seg_length         <- 10
n_seg_per_transect <- floor(X_km_max / seg_length)

transect_Y <- seq(Y_km_max / (2 * n_transects),
                  Y_km_max - Y_km_max / (2 * n_transects),
                  length.out = n_transects)

segments <- expand.grid(
  transect_id = seq_len(n_transects),
  seg_idx     = seq_len(n_seg_per_transect)
) |>
  dplyr::mutate(
    X_mid    = (seg_idx - 0.5) * seg_length,
    Y_mid    = transect_Y[transect_id],
    Z_bathy  = bathy_mean_fn(X_mid, Y_mid),
    seg_id   = dplyr::row_number(),
    seg_l    = seg_length
  )

n_seg_total <- nrow(segments)
total_effort_km <- sum(segments$seg_l)

cat(sprintf("Survey: %d transects x %d segments = %d segments\n",
            n_transects, n_seg_per_transect, n_seg_total))
cat(sprintf("  total effort = %.0f km, surveyed area = %.0f km^2\n",
            total_effort_km, total_effort_km * 2 * w))
cat(sprintf("  truncation distance w = %.2f km\n\n", w))

# Coords at segment midpoints, normalised to [-1, 1]
coords_raw <- as.matrix(segments[, c("X_mid", "Y_mid", "Z_bathy")])
coords_norm <- sweep(coords_raw, 2, coord_centre, "-")
coords_norm <- sweep(coords_norm, 2, coord_scale,  "/")
stopifnot(all(is.finite(coords_norm)))
stopifnot(max(abs(coords_norm)) < 1.5)

cat(sprintf("Coord normalisation: centre=[%.0f,%.0f,%.0f], scale=[%.0f,%.0f,%.0f]\n",
            coord_centre[1], coord_centre[2], coord_centre[3],
            coord_scale[1], coord_scale[2], coord_scale[3]))
cat(sprintf("  X_norm: [%.3f, %.3f]\n",
            min(coords_norm[, 1]), max(coords_norm[, 1])))
cat(sprintf("  Y_norm: [%.3f, %.3f]\n",
            min(coords_norm[, 2]), max(coords_norm[, 2])))
cat(sprintf("  Z_norm: [%.3f, %.3f]\n",
            min(coords_norm[, 3]), max(coords_norm[, 3])))

# -----------------------------------------------------------------------------
# 5. Per-species simulation + fit
# -----------------------------------------------------------------------------

stan_file <- "distance/distance_v4.1a.stan"
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

iter_warmup   <- as.integer(Sys.getenv("ITER_WARMUP",   "1000"))
iter_sampling <- as.integer(Sys.getenv("ITER_SAMPLING", "1000"))
max_treedepth <- as.integer(Sys.getenv("MAX_TREEDEPTH", "12"))

cat(sprintf("\nStan settings: warmup=%d, sampling=%d, max_treedepth=%d\n",
            iter_warmup, iter_sampling, max_treedepth))

results <- list()

for (sp_key in names(species_params)) {
  p <- species_params[[sp_key]]
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("=== %s\n", p$common_name))
  cat(sprintf("======================================================\n"))

  # ---- 5a. Draw GP at segment midpoints (simulation truth) -----------------
  K     <- aniso_cov(coords_raw, p$sigma, p$lx, p$ly, p$lz)
  f_seg <- as.vector(MASS::mvrnorm(1, mu = rep(0, n_seg_total), Sigma = K))

  # In this spatial model, mu_sp is the log GROUP density intercept.
  # Truth at segments: log_lambda_groups = log(exp(mu)/mean_gs) + f_seg
  mu_groups_truth   <- p$mu - log(p$mean_group_size)
  log_lambda_groups <- mu_groups_truth + f_seg
  lambda_groups     <- exp(log_lambda_groups)
  lambda_animals    <- lambda_groups * p$mean_group_size

  cat(sprintf("\n  spatial truth (mean across %d segments):\n", n_seg_total))
  cat(sprintf("    mu_sp_truth     = %.3f (log groups/km^2)\n", mu_groups_truth))
  cat(sprintf("    lambda_groups   = %.5f (mean), %.5f (median)\n",
              mean(lambda_groups), median(lambda_groups)))
  cat(sprintf("    lambda_animals  = %.4f (mean), %.4f (median)\n",
              mean(lambda_animals), median(lambda_animals)))
  cat(sprintf("    f_seg range     = [%.2f, %.2f]\n",
              min(f_seg), max(f_seg)))

  # ---- 5b. Simulate group sightings (verbatim from v4.1) -------------------
  expected_groups   <- lambda_groups * 2 * w * seg_length
  n_groups_in_strip <- rpois(n_seg_total, expected_groups)

  detect_rows <- vector("list", n_seg_total)
  for (i in seq_len(n_seg_total)) {
    n_strip <- n_groups_in_strip[i]
    if (n_strip == 0L) next
    x_true <- runif(n_strip, 0, w)
    if (p$sim_group_dist == "empirical") {
      s_strip <- sample(p$group_size_pool, n_strip, replace = TRUE)
    } else if (p$sim_group_dist == "constant") {
      s_strip <- rep(p$mean_group_size, n_strip)
    } else if (p$sim_group_dist == "nb") {
      s_strip <- rnbinom_zt(n_strip,
                            mu   = p$group_size_nb$mu,
                            size = p$group_size_nb$size)
    }
    if (p$use_size_covar == 1L) {
      sigma_g <- p$sigma_det * exp(p$beta_size_truth * (s_strip - p$s_centre))
    } else {
      sigma_g <- rep(p$sigma_det, n_strip)
    }
    p_det <- exp(-x_true^2 / (2 * sigma_g^2))
    keep  <- runif(n_strip) < p_det
    if (any(keep)) {
      detect_rows[[i]] <- data.frame(
        seg_id   = segments$seg_id[i],
        distance = x_true[keep],
        size     = s_strip[keep]
      )
    }
  }
  obs <- do.call(rbind, detect_rows)
  if (is.null(obs)) {
    obs <- data.frame(seg_id = integer(), distance = numeric(), size = integer())
  }

  seg_count <- tabulate(obs$seg_id, nbins = n_seg_total)

  cat(sprintf("\n  groups in strip (truth): %d\n", sum(n_groups_in_strip)))
  cat(sprintf("  groups detected         : %d (%.1f%% of strip)\n",
              nrow(obs), 100 * nrow(obs) / max(sum(n_groups_in_strip), 1L)))
  cat(sprintf("  segments with >= 1 detection: %d / %d\n",
              sum(seg_count > 0), n_seg_total))

  if (nrow(obs) < 10L) {
    warning(sprintf("Too few detections (%d) for %s - skipping fit.",
                    nrow(obs), p$common_name))
    next
  }

  # ---- 5c. Stan data list ---------------------------------------------------
  stan_data <- list(

    # Detection-level data
    n                       = nrow(obs),
    x                       = obs$distance,
    s                       = as.integer(obs$size),
    det_seg                 = as.integer(obs$seg_id),
    w                       = w,

    # Detection function priors
    log_sigma_prior_mean    = p$log_sigma_prior_mean,
    log_sigma_prior_sd      = p$log_sigma_prior_sd,
    use_size_covar          = p$use_size_covar,
    s_centre                = p$s_centre,
    beta_size_prior_mean    = p$beta_size_prior_mean,
    beta_size_prior_sd      = p$beta_size_prior_sd,
    S_max                   = p$S_max,

    # Group-size sub-model
    group_size_dist         = p$model_group_dist,
    mu_s_prior_shape        = p$mu_s_prior_shape,
    mu_s_prior_rate         = p$mu_s_prior_rate,
    phi_s_prior_shape       = p$phi_s_prior_shape,
    phi_s_prior_rate        = p$phi_s_prior_rate,
    mu_log_prior_mean       = p$mu_log_prior_mean,
    mu_log_prior_sd         = p$mu_log_prior_sd,
    sigma_log_prior_shape   = p$sigma_log_prior_shape,
    sigma_log_prior_rate    = p$sigma_log_prior_rate,

    # Per-segment encounter-rate data
    n_seg                   = n_seg_total,
    seg_count               = as.integer(seg_count),
    seg_l                   = rep(seg_length, n_seg_total),

    # Spatial / HSGP block
    S                       = 1L,
    M                       = prod(HSGP_M),
    coords                  = coords_norm,
    coord_scale             = coord_scale,
    L_hsgp                  = HSGP_C,
    m_hsgp                  = HSGP_M,
    prior_mu_sp_mu          = p$prior_mu_sp_mu,
    prior_mu_sp_sig         = p$prior_mu_sp_sig,
    prior_gp_sigma_shape    = p$prior_gp_sigma_shape,
    prior_gp_sigma_rate     = p$prior_gp_sigma_rate,
    prior_gp_lx_mu          = p$prior_gp_lx_mu,
    prior_gp_lx_sig         = p$prior_gp_lx_sig,
    prior_gp_ly_mu          = p$prior_gp_ly_mu,
    prior_gp_ly_sig         = p$prior_gp_ly_sig,
    prior_gp_lz_mu          = p$prior_gp_lz_mu,
    prior_gp_lz_sig         = p$prior_gp_lz_sig
  )

  # ---- 5d. Init function ----------------------------------------------------
  init_fn <- function() {
    list(
      mu_sp        = array(p$prior_mu_sp_mu + rnorm(1, 0, 0.3), dim = 1),
      gp_sigma     = array(p$sigma + abs(rnorm(1, 0, 0.2)), dim = 1),
      gp_l         = matrix(c(p$prior_gp_lx_mu, p$prior_gp_ly_mu,
                              p$prior_gp_lz_mu),
                            nrow = 1) +
                     matrix(rnorm(3, 0, 5), nrow = 1),
      z_beta       = matrix(rnorm(prod(HSGP_M) * 1, 0, 0.1),
                            nrow = 1),
      log_sigma    = p$log_sigma_prior_mean + rnorm(1, 0, 0.3),
      beta_size    = rnorm(1, 0, 0.01),
      mu_s         = max(p$mu_s_prior_shape / p$mu_s_prior_rate
                         + rnorm(1, 0, 0.5 * p$mu_s_prior_shape /
                                          p$mu_s_prior_rate),
                         0.5),
      phi_s        = runif(1, 0.5, 2),
      mu_log_s     = p$mu_log_prior_mean + rnorm(1, 0, 0.2),
      sigma_log_s  = max(p$sigma_log_prior_shape / p$sigma_log_prior_rate
                         + rnorm(1, 0, 0.2),
                         0.1)
    )
  }

  # ---- 5e. Stan fit ---------------------------------------------------------
  cat("\n  Fitting distance_v4.1a.stan (HSGP spatial)...\n")
  fit <- mod$sample(
    data              = stan_data,
    chains            = 4,
    parallel_chains   = 4,
    threads_per_chain = 1,
    iter_warmup       = iter_warmup,
    iter_sampling     = iter_sampling,
    max_treedepth     = max_treedepth,
    seed              = 42,
    init              = init_fn,
    refresh           = 200,
    show_messages     = FALSE
  )
  invisible(fit$cmdstan_diagnose())

  draws <- fit$draws(format = "df")

  # ---- 5f. Spatial parameter recovery --------------------------------------
  # Per-segment posterior summary of lambda_groups (transformed parameter)
  lambda_groups_arr <- fit$draws("lambda_groups", format = "matrix")
  f_seg_arr         <- fit$draws("f_seg",         format = "matrix")
  seg_count_rep_arr <- fit$draws("seg_count_rep", format = "matrix")

  lambda_groups_med <- apply(lambda_groups_arr, 2, median)
  lambda_groups_lo  <- apply(lambda_groups_arr, 2, quantile, 0.025)
  lambda_groups_hi  <- apply(lambda_groups_arr, 2, quantile, 0.975)

  f_seg_med <- apply(f_seg_arr, 2, median)
  f_seg_lo  <- apply(f_seg_arr, 2, quantile, 0.025)
  f_seg_hi  <- apply(f_seg_arr, 2, quantile, 0.975)

  seg_count_rep_med <- apply(seg_count_rep_arr, 2, median)

  # Coverage on per-segment lambda_groups
  cover_lambda <- mean(
    lambda_groups >= lambda_groups_lo &
    lambda_groups <= lambda_groups_hi
  )
  cover_f <- mean(f_seg >= f_seg_lo & f_seg <= f_seg_hi)
  cor_lambda <- cor(log(lambda_groups), log(lambda_groups_med))
  cor_f      <- cor(f_seg, f_seg_med)

  cat("\n  Spatial recovery (per-segment):\n")
  cat(sprintf("    log(lambda_groups): truth-vs-posterior cor = %.3f, 95%% CI cover = %.1f%%\n",
              cor_lambda, 100 * cover_lambda))
  cat(sprintf("    f_seg            : truth-vs-posterior cor = %.3f, 95%% CI cover = %.1f%%\n",
              cor_f,      100 * cover_f))

  # GP hyperparameter recovery
  gp_recovery <- tibble(
    param  = c("mu_sp", "gp_sigma", "lx", "ly", "lz"),
    truth  = c(mu_groups_truth,
               p$sigma,
               p$lx, p$ly, p$lz),
    median = c(median(draws$`mu_sp[1]`),
               median(draws$`gp_sigma[1]`),
               median(draws$`gp_l[1,1]`),
               median(draws$`gp_l[1,2]`),
               median(draws$`gp_l[1,3]`)),
    q025   = c(quantile(draws$`mu_sp[1]`,    0.025),
               quantile(draws$`gp_sigma[1]`, 0.025),
               quantile(draws$`gp_l[1,1]`,   0.025),
               quantile(draws$`gp_l[1,2]`,   0.025),
               quantile(draws$`gp_l[1,3]`,   0.025)),
    q975   = c(quantile(draws$`mu_sp[1]`,    0.975),
               quantile(draws$`gp_sigma[1]`, 0.975),
               quantile(draws$`gp_l[1,1]`,   0.975),
               quantile(draws$`gp_l[1,2]`,   0.975),
               quantile(draws$`gp_l[1,3]`,   0.975))
  )
  cat("\n  GP hyperparameter recovery:\n")
  print(knitr::kable(gp_recovery, digits = 3, format = "simple"))

  # Detection-side recovery (same structure as v4.1)
  if (p$model_group_dist == 0L) {
    p_zero_post     <- (draws$phi_s / (draws$mu_s + draws$phi_s))^draws$phi_s
    mean_s_obs_draw <- draws$mu_s / (1 - p_zero_post)
  } else {
    mean_s_obs_draw <- exp(draws$mu_log_s + 0.5 * draws$sigma_log_s^2)
  }

  det_recovery <- tibble(
    param = c("sigma (km)",
              "beta_size",
              "E[s|s>0]"),
    truth = c(p$sigma_det,
              if (p$use_size_covar == 1L) p$beta_size_truth else NA_real_,
              p$mean_group_size),
    median = c(median(exp(draws$log_sigma)),
               median(draws$beta_size),
               median(mean_s_obs_draw)),
    q025   = c(quantile(exp(draws$log_sigma), 0.025),
               quantile(draws$beta_size, 0.025),
               quantile(mean_s_obs_draw, 0.025)),
    q975   = c(quantile(exp(draws$log_sigma), 0.975),
               quantile(draws$beta_size, 0.975),
               quantile(mean_s_obs_draw, 0.975))
  )
  cat("\n  Detection-side recovery:\n")
  print(knitr::kable(det_recovery, digits = 4, format = "simple"))

  # ---- 5g. Plots ------------------------------------------------------------
  sp_colour <- if (sp_key == "humpback") "#3B528B" else "#5DC863"

  # (A) Truth lambda surface at segment midpoints
  p_truth <- ggplot(
      segments |> dplyr::mutate(lambda = lambda_groups),
      aes(X_mid, Y_mid, colour = lambda)) +
    geom_point(size = 1.0, alpha = 0.85) +
    scale_colour_viridis_c(name = expression(lambda[truth]),
                           trans = "log10") +
    coord_fixed(xlim = c(0, X_km_max), ylim = c(0, Y_km_max), expand = FALSE) +
    labs(title    = "Simulated TRUTH",
         subtitle = "Per-segment lambda_groups (groups/km^2)",
         x = "Easting (km)", y = "Northing (km)") +
    theme_bw(base_size = 10)

  # (B) Posterior median lambda surface
  p_post <- ggplot(
      segments |> dplyr::mutate(lambda = lambda_groups_med),
      aes(X_mid, Y_mid, colour = lambda)) +
    geom_point(size = 1.0, alpha = 0.85) +
    scale_colour_viridis_c(name = expression(lambda[post]),
                           trans = "log10") +
    coord_fixed(xlim = c(0, X_km_max), ylim = c(0, Y_km_max), expand = FALSE) +
    labs(title    = "POSTERIOR median",
         subtitle = "Per-segment lambda_groups (groups/km^2)",
         x = "Easting (km)", y = NULL) +
    theme_bw(base_size = 10)

  # (C) Per-segment truth-vs-posterior scatter for log(lambda_groups)
  p_lambda_scatter <- ggplot(
      tibble(truth = log(lambda_groups),
             post  = log(lambda_groups_med)),
      aes(truth, post)) +
    geom_abline(slope = 1, intercept = 0,
                colour = "grey50", linetype = "dashed") +
    geom_point(colour = sp_colour, size = 1.4, alpha = 0.6) +
    labs(title    = "PPC: per-segment lambda_groups",
         subtitle = sprintf("cor(log truth, log posterior median) = %.3f, 95%% CI cover = %.1f%%",
                            cor_lambda, 100 * cover_lambda),
         x = expression(log(lambda[truth])),
         y = expression(log(lambda[posterior]))) +
    theme_bw(base_size = 10)

  # (D) Per-segment truth-vs-posterior scatter for f_seg
  p_f_scatter <- ggplot(
      tibble(truth = f_seg,
             post  = f_seg_med),
      aes(truth, post)) +
    geom_abline(slope = 1, intercept = 0,
                colour = "grey50", linetype = "dashed") +
    geom_point(colour = sp_colour, size = 1.4, alpha = 0.6) +
    labs(title    = "PPC: spatial component f_seg",
         subtitle = sprintf("cor(truth, posterior median) = %.3f, 95%% CI cover = %.1f%%",
                            cor_f, 100 * cover_f),
         x = "f_seg (truth)",
         y = "f_seg (posterior)") +
    theme_bw(base_size = 10)

  # (E) PPC: posterior-predictive segment counts vs observed
  ppc_df <- tibble(
    seg_id        = segments$seg_id,
    obs_count     = seg_count,
    rep_med       = seg_count_rep_med,
    rep_lo        = apply(seg_count_rep_arr, 2, quantile, 0.025),
    rep_hi        = apply(seg_count_rep_arr, 2, quantile, 0.975)
  )
  p_count_ppc <- ggplot(ppc_df, aes(obs_count, rep_med)) +
    geom_abline(slope = 1, intercept = 0,
                colour = "grey50", linetype = "dashed") +
    geom_errorbar(aes(ymin = rep_lo, ymax = rep_hi),
                  colour = sp_colour, alpha = 0.25, width = 0) +
    geom_point(colour = sp_colour, size = 1.4, alpha = 0.7) +
    labs(title    = "PPC: posterior-predictive segment counts",
         subtitle = sprintf("Observed total = %d; posterior-pred mean = %.1f",
                            sum(seg_count), mean(rowSums(seg_count_rep_arr))),
         x = "Observed seg_count",
         y = "Posterior-predictive seg_count_rep (median, 95% CI)") +
    theme_bw(base_size = 10)

  # (F) GP hyperparameter recovery scatter
  gp_plot_df <- gp_recovery |>
    dplyr::filter(!is.na(truth)) |>
    dplyr::mutate(param = factor(param, levels = param))
  p_gp_recov <- ggplot(gp_plot_df,
                       aes(x = param, y = median,
                           ymin = q025, ymax = q975)) +
    geom_errorbar(width = 0.15, colour = sp_colour, linewidth = 0.8) +
    geom_point(colour = sp_colour, size = 3) +
    geom_point(aes(y = truth), colour = "red", shape = 4,
               size = 3, stroke = 1) +
    facet_wrap(~param, scales = "free", nrow = 1) +
    labs(title    = "GP hyperparameter recovery",
         subtitle = "Posterior median +/- 95% CI; red x = truth",
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())

  combined <- (p_truth | p_post) /
              (p_lambda_scatter | p_f_scatter) /
              (p_count_ppc | p_gp_recov) +
              plot_annotation(
                title    = sprintf("Distance-sampling spatial HSGP (v4.1a) -- %s",
                                   p$common_name),
                subtitle = sprintf("HN(sigma=%.2f km, use_size_covar=%d), w=%.2f km. %d transects x %d segments x 10 km. %d detections in %d segments. M_HSGP = (%d, %d, %d).",
                                   p$sigma_det, p$use_size_covar, w,
                                   n_transects, n_seg_per_transect,
                                   nrow(obs), sum(seg_count > 0),
                                   HSGP_M[1], HSGP_M[2], HSGP_M[3]),
                theme = theme(plot.title    = element_text(face = "bold",
                                                           size = 13),
                              plot.subtitle = element_text(size = 9,
                                                           colour = "grey40"))
              )

  out_pdf <- file.path(OUTPUT_DIR,
                       sprintf("distance_v4.1a_%s.pdf", sp_key))
  grDevices::cairo_pdf(out_pdf, width = 12, height = 13)
  print(combined)
  invisible(dev.off())
  cat(sprintf("\n  Saved %s\n", out_pdf))

  # Save fit + sim metadata
  saveRDS(list(
    sp_key            = sp_key,
    sp_params         = p,
    w                 = w,
    transect_Y        = transect_Y,
    seg_length        = seg_length,
    segments          = segments,
    coords_norm       = coords_norm,
    f_seg             = f_seg,
    lambda_groups     = lambda_groups,
    lambda_animals    = lambda_animals,
    n_groups_in_strip = n_groups_in_strip,
    obs               = obs,
    seg_count         = seg_count,
    stan_data         = stan_data,
    fit               = fit,
    gp_recovery       = gp_recovery,
    det_recovery      = det_recovery,
    cover_lambda      = cover_lambda,
    cover_f           = cover_f
  ), file.path(OUTPUT_DIR, sprintf("distance_v4.1a_%s.rds", sp_key)))

  results[[sp_key]] <- bind_rows(
    gp_recovery  |> dplyr::mutate(species = p$common_name, group = "GP"),
    det_recovery |> dplyr::mutate(species = p$common_name, group = "detection")
  )
}

# -----------------------------------------------------------------------------
# 6. Combined recovery table
# -----------------------------------------------------------------------------
if (length(results) > 0L) {
  cat("\n\n======================================================\n")
  cat("=== Combined recovery table (GP + detection)\n")
  cat("======================================================\n\n")

  combined_summary <- do.call(rbind, results) |>
    dplyr::select(species, group, param, truth, median, q025, q975)
  print(knitr::kable(combined_summary, digits = 4, format = "simple"))
  write_csv(combined_summary,
            file.path(OUTPUT_DIR, "distance_v4.1a_recovery.csv"))
}
