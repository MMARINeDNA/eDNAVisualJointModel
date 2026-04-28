# =============================================================================
# 00_distance_v4.1.R
#
# Test simulation pipeline for line-transect distance sampling, mirroring
# the v4.1 spatial-sim conventions for the underlying density field.
#
# Per species (humpback whale and Pacific white-sided dolphin), this
# script:
#
#   1. Simulates a zero-mean GP density surface (animals/km^2) using the
#      v4.1 domain, bathymetry, and species-specific GP hyperparameters.
#   2. Converts the animal-density surface to a group-density surface:
#         lambda_groups = lambda_animals / mean_group_size,
#      with mean group sizes 2 (humpback) and 50 (PWSD).
#   3. Lays down systematic parallel north-south transects, splits them
#      into 10 km segments, and simulates group sightings from a
#      half-normal detection function with truncation distance w.
#   4. Assigns a constant group size (the species mean) to every
#      detected group - a placeholder until we add a group-size model.
#   5. Fits distance/distance_hn_dens.stan to the simulated
#      observations (one fit per species).
#   6. Writes diagnostic plots and a parameter-recovery table to
#      outputs/distance_v4.1/.
#
# The simulated data are drawn from a SPATIAL field but fit by a
# NON-SPATIAL distance-sampling model (lambda_s as a single scalar).
# So the recovery target for lambda_s is the spatial mean of the
# species' group-density surface, not any local value.
#
# Run from the project root:
#   Rscript distance/00_distance_v4.1.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(cmdstanr)
  library(posterior)
  library(patchwork)
  library(viridis)
})

set.seed(42)

OUTPUT_DIR <- "outputs/distance_v4.1"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 0. Helpers
# -----------------------------------------------------------------------------

# Standard error function in terms of the normal CDF (avoids extra package).
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# Anisotropic squared-exponential covariance with diagonal jitter.
aniso_cov <- function(coords, sigma, lx, ly, lz, jitter = 1e-6) {
  d2 <- outer(coords[, 1], coords[, 1], function(a, b) ((a - b) / lx)^2) +
        outer(coords[, 2], coords[, 2], function(a, b) ((a - b) / ly)^2) +
        outer(coords[, 3], coords[, 3], function(a, b) ((a - b) / lz)^2)
  sigma * sigma * exp(-0.5 * d2) + diag(jitter, nrow(d2))
}

# Anisotropic SE cross-covariance, for kriging onto a separate grid.
aniso_cross_cov <- function(c1, c2, sigma, lx, ly, lz) {
  d2 <- outer(c1[, 1], c2[, 1], function(a, b) ((a - b) / lx)^2) +
        outer(c1[, 2], c2[, 2], function(a, b) ((a - b) / ly)^2) +
        outer(c1[, 3], c2[, 3], function(a, b) ((a - b) / lz)^2)
  sigma * sigma * exp(-0.5 * d2)
}

# -----------------------------------------------------------------------------
# 1. Domain (matches v4.1)
# -----------------------------------------------------------------------------

X_km_max     <- 500          # E-W extent (km)
Y_km_max     <- 1270         # N-S extent (km)
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

# -----------------------------------------------------------------------------
# 2. Species (humpback + PWSD)
# -----------------------------------------------------------------------------

species_params <- list(
  humpback = list(
    common_name      = "Humpback whale",
    sigma            = 1.0,
    lx               = 50,        # km
    ly               = 300,       # km
    lz               = 100,       # m
    mu               = log(0.004),# log animals/km^2
    mean_group_size  = 2L
  ),
  pwsd = list(
    common_name      = "Pacific white-sided dolphin",
    sigma            = 1.3,
    lx               = 40,
    ly               = 300,
    lz               = 300,
    mu               = log(0.032),
    mean_group_size  = 50L
  )
)

# -----------------------------------------------------------------------------
# 3. Detection function (shared across species)
# -----------------------------------------------------------------------------

sigma_det <- 2.0    # half-normal scale (km)
w         <- 5.0    # truncation distance (km)
esw       <- sigma_det * sqrt(pi / 2) * erf(w / (sigma_det * sqrt(2)))
p_avg     <- esw / w
cat(sprintf("Detection function: HN(sigma=%.2f km), truncation w=%.2f km\n",
            sigma_det, w))
cat(sprintf("  effective strip half-width ESW = %.3f km\n", esw))
cat(sprintf("  mean detection prob in strip   = %.3f\n", p_avg))

# -----------------------------------------------------------------------------
# 4. Survey design: systematic parallel E-W transects, 10 km segments
#
# Transects run E-W (variable X, fixed Y). With Y_km_max = 1270 km and
# n_transects, transect spacing in Y is Y_km_max / n_transects, offset
# from each edge by half-spacing. Each transect spans X from 0 to
# X_km_max = 500 km, split into 10-km segments.
# -----------------------------------------------------------------------------

n_transects        <- 25                         # bumped to keep coverage
seg_length         <- 10                         # km
n_seg_per_transect <- floor(X_km_max / seg_length)

# Transect Y-positions, offset by half-spacing from each edge
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

cat(sprintf("\nSurvey: %d transects x %d segments = %d segments\n",
            n_transects, n_seg_per_transect, n_seg_total))
cat(sprintf("  total effort = %.0f km, surveyed area = %.0f km^2\n",
            total_effort_km, total_effort_km * 2 * w))

# -----------------------------------------------------------------------------
# 5. Per-species simulation + fit
# -----------------------------------------------------------------------------

stan_file <- "distance/distance_hn_dens.stan"
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

results <- list()

for (sp_key in names(species_params)) {
  p <- species_params[[sp_key]]
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("=== %s\n", p$common_name))
  cat(sprintf("======================================================\n"))

  # ------------- 5a. Draw GP at segment midpoints -----------------------------
  coords <- as.matrix(segments[, c("X_mid", "Y_mid", "Z_bathy")])
  K      <- aniso_cov(coords, p$sigma, p$lx, p$ly, p$lz)
  f_seg  <- as.vector(MASS::mvrnorm(1, mu = rep(0, n_seg_total), Sigma = K))

  lambda_animals <- exp(p$mu + f_seg)                          # animals/km^2
  lambda_groups  <- lambda_animals / p$mean_group_size         # groups/km^2

  cat(sprintf("\n  spatial truth (mean across %d segments):\n", n_seg_total))
  cat(sprintf("    lambda_animals = %.4f animals/km^2\n", mean(lambda_animals)))
  cat(sprintf("    lambda_groups  = %.4f groups/km^2\n",  mean(lambda_groups)))

  # ------------- 5b. Simulate group sightings ---------------------------------

  # Number of groups in each segment's strip (Poisson)
  expected_groups   <- lambda_groups * 2 * w * seg_length
  n_groups_in_strip <- rpois(n_seg_total, expected_groups)

  # For each group: uniform distance in [0, w], thin by HN detection
  detect_rows <- vector("list", n_seg_total)
  for (i in seq_len(n_seg_total)) {
    n_strip <- n_groups_in_strip[i]
    if (n_strip == 0L) next
    x_true <- runif(n_strip, 0, w)
    p_det  <- exp(-x_true^2 / (2 * sigma_det^2))
    keep   <- runif(n_strip) < p_det
    if (any(keep)) {
      detect_rows[[i]] <- data.frame(
        seg_id   = segments$seg_id[i],
        distance = x_true[keep],
        size     = p$mean_group_size
      )
    }
  }
  obs <- do.call(rbind, detect_rows)
  if (is.null(obs)) {
    obs <- data.frame(seg_id = integer(), distance = numeric(), size = integer())
  }

  seg_count <- tabulate(obs$seg_id, nbins = n_seg_total)

  cat(sprintf("\n  groups in strip (truth) : %d\n", sum(n_groups_in_strip)))
  cat(sprintf("  groups detected         : %d (%.1f%% of strip)\n",
              nrow(obs), 100 * nrow(obs) / max(sum(n_groups_in_strip), 1L)))
  cat(sprintf("  segments with >= 1 detection: %d / %d\n",
              sum(seg_count > 0), n_seg_total))

  # ------------- 5c. Stan data ------------------------------------------------

  if (nrow(obs) < 10L) {
    warning(sprintf("Too few detections (%d) for %s - skipping fit.",
                    nrow(obs), p$common_name))
    next
  }

  stan_data <- list(
    n                    = nrow(obs),
    x                    = obs$distance,
    s                    = as.integer(obs$size),
    w                    = w,
    log_sigma_prior_mean = log(w / 2),
    log_sigma_prior_sd   = 0.5,
    n_seg                = n_seg_total,
    seg_count            = as.integer(seg_count),
    seg_l                = rep(seg_length, n_seg_total)
  )

  # ------------- 5d. Stan fit -------------------------------------------------

  init_fn <- function() {
    lam_init <- max(mean(stan_data$seg_count /
                         (stan_data$seg_l * 2 * w)), 1e-6)
    list(
      log_sigma    = log(runif(1, w / 100, w * 2)),
      mu_s         = runif(1, 1, 100),
      phi_s        = runif(1, 0.1, 1),
      log_lambda_s = log(runif(1, lam_init * 0.1, lam_init * 1.9 + 1e-6))
    )
  }

  cat("\n  Fitting distance_hn_dens.stan...\n")
  fit <- mod$sample(
    data              = stan_data,
    chains            = 4,
    parallel_chains   = 4,
    threads_per_chain = 1,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    seed              = 42,
    init              = init_fn,
    refresh           = 500,
    show_messages     = FALSE
  )

  diag <- fit$cmdstan_diagnose()

  # ------------- 5e. Extract draws + summarise --------------------------------

  draws <- fit$draws(format = "df")
  D_animals_draw <- draws$lambda_s * draws$mu_s

  param_summary <- tibble(
    param      = c("sigma (km)", "mu_s (mean group size)",
                   "lambda_s (groups/km^2)", "D (animals/km^2)"),
    truth      = c(sigma_det, p$mean_group_size,
                   mean(lambda_groups), mean(lambda_animals)),
    median     = c(median(draws$sigma), median(draws$mu_s),
                   median(draws$lambda_s), median(D_animals_draw)),
    q025       = c(quantile(draws$sigma, 0.025), quantile(draws$mu_s, 0.025),
                   quantile(draws$lambda_s, 0.025),
                   quantile(D_animals_draw, 0.025)),
    q975       = c(quantile(draws$sigma, 0.975), quantile(draws$mu_s, 0.975),
                   quantile(draws$lambda_s, 0.975),
                   quantile(D_animals_draw, 0.975))
  )

  cat("\n  Parameter recovery (truth = spatial mean for lambda / D):\n")
  print(knitr::kable(param_summary, digits = 4, format = "simple"))

  # ------------- 5f. Plots ----------------------------------------------------

  sp_colour <- if (sp_key == "humpback") "#3B528B" else "#5DC863"

  # (A) Density surface scatter at segment midpoints
  p_density <- ggplot(segments |> mutate(lambda_animals = lambda_animals,
                                         lambda_groups  = lambda_groups),
                      aes(X_mid, Y_mid, colour = lambda_animals)) +
    geom_point(size = 1.0, alpha = 0.85) +
    scale_colour_viridis_c(name = expression(lambda ~ "(animals/km"^2 * ")"),
                           trans = "log10") +
    coord_fixed(xlim = c(0, X_km_max), ylim = c(0, Y_km_max), expand = FALSE) +
    labs(title    = sprintf("%s: simulated density at segment midpoints",
                            p$common_name),
         subtitle = "Each point = one 10-km transect segment, coloured by truth lambda",
         x = "Easting (km)", y = "Northing (km)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))

  # (B) Detection-distance histogram + fitted HN
  sigma_post <- median(draws$sigma)
  esw_post   <- sigma_post * sqrt(pi / 2) * erf(w / (sigma_post * sqrt(2)))

  p_hist <- ggplot(obs, aes(x = distance)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "lightblue", colour = "white") +
    stat_function(fun = function(x) {
      exp(-x^2 / (2 * sigma_post^2)) / esw_post
    }, colour = "red", linewidth = 1) +
    stat_function(fun = function(x) {
      exp(-x^2 / (2 * sigma_det^2)) / esw
    }, colour = "darkgrey", linetype = "dashed", linewidth = 0.8) +
    labs(title    = sprintf("%s: detection function recovery", p$common_name),
         subtitle = sprintf("Red = posterior median HN (sigma=%.2f); grey dashed = truth (sigma=%.2f)",
                            sigma_post, sigma_det),
         x = "Detection distance (km)", y = "Density") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))

  # (C) Posterior intervals vs truth, log-transformed where useful
  param_summary_long <- param_summary |>
    dplyr::mutate(param = factor(param, levels = param))

  p_recovery <- ggplot(param_summary_long,
                       aes(x = param,
                           y = median,
                           ymin = q025, ymax = q975)) +
    geom_errorbar(width = 0.15, colour = sp_colour, linewidth = 0.8) +
    geom_point(colour = sp_colour, size = 3) +
    geom_point(aes(y = truth), colour = "red", shape = 4, size = 3, stroke = 1) +
    facet_wrap(~param, scales = "free", nrow = 1) +
    labs(title    = sprintf("%s: parameter recovery", p$common_name),
         subtitle = "Posterior median +/- 95% CI (coloured); red x = truth",
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title  = element_text(face = "bold", size = 11),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text  = element_text(size = 9))

  # (D) Posterior of D (animals/km^2) overlay
  truth_D <- mean(lambda_animals)
  p_D <- ggplot(tibble(D = D_animals_draw), aes(x = D)) +
    geom_density(fill = sp_colour, colour = NA, alpha = 0.4) +
    geom_density(colour = sp_colour, linewidth = 0.9) +
    geom_vline(xintercept = truth_D, colour = "red",
               linetype = "dashed", linewidth = 0.7) +
    labs(title    = sprintf("%s: posterior of D (animals/km^2)",
                            p$common_name),
         subtitle = sprintf("Red dashed = simulated truth (spatial mean lambda) = %.4f",
                            truth_D),
         x = expression(D ~ "(animals/km"^2 * ")"),
         y = "posterior density") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))

  combined <- (p_density / (p_hist | p_D)) /
              p_recovery +
              plot_layout(heights = c(2, 1.4, 1)) +
              plot_annotation(
                title    = sprintf("Distance sampling pipeline (v4.1) -- %s",
                                   p$common_name),
                subtitle = sprintf("HN(sigma=%.2f km), w=%.2f km. %d transects x %d segments x 10 km. %d detections.",
                                   sigma_det, w, n_transects, n_seg_per_transect, nrow(obs)),
                theme = theme(plot.title    = element_text(face = "bold", size = 13),
                              plot.subtitle = element_text(size = 10, colour = "grey40"))
              )

  out_pdf <- file.path(OUTPUT_DIR, sprintf("distance_v4.1_%s.pdf", sp_key))
  grDevices::cairo_pdf(out_pdf, width = 11, height = 11)
  print(combined)
  invisible(dev.off())
  cat(sprintf("  Saved %s\n", out_pdf))

  # Save fit + sim metadata
  saveRDS(list(
    sp_key            = sp_key,
    sp_params         = p,
    sigma_det         = sigma_det,
    w                 = w,
    transect_Y        = transect_Y,
    seg_length        = seg_length,
    segments          = segments,
    f_seg             = f_seg,
    lambda_animals    = lambda_animals,
    lambda_groups     = lambda_groups,
    n_groups_in_strip = n_groups_in_strip,
    obs               = obs,
    seg_count         = seg_count,
    stan_data         = stan_data,
    fit               = fit,
    param_summary     = param_summary
  ), file.path(OUTPUT_DIR, sprintf("distance_v4.1_%s.rds", sp_key)))

  results[[sp_key]] <- param_summary |> dplyr::mutate(species = p$common_name)
}

# -----------------------------------------------------------------------------
# 6. Combined summary
# -----------------------------------------------------------------------------

if (length(results) > 0) {
  cat("\n\n======================================================\n")
  cat("=== Combined recovery table\n")
  cat("======================================================\n\n")

  combined_summary <- do.call(rbind, results) |>
    dplyr::select(species, param, truth, median, q025, q975)
  print(knitr::kable(combined_summary, digits = 4, format = "simple"))
  write_csv(combined_summary,
            file.path(OUTPUT_DIR, "distance_v4.1_recovery.csv"))
}

cat(sprintf("\n=== Done. Outputs written to %s ===\n", OUTPUT_DIR))
for (f in list.files(OUTPUT_DIR)) cat(sprintf("  %s\n", f))
