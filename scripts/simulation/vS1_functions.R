# =============================================================================
# vS1_functions.R
#
# Callable versions of the vS1 simulation (01) and Stan-data formatting (03)
# stages, parameterised by RNG seed and HSGP basis count. Used by the basis
# sweep (06_basis_sweep_vS1.r); the canonical single-run scripts 01/03 are
# left untouched.
#
# Both functions mirror the committed 01 / 03 logic (including the
# junk-constant-within-sample fix). simulate_whale_edna_vS1() differs from 01
# in one intended way: it seeds the whole draw from a single `seed` argument
# (01 uses several fixed internal seeds to pin one reference dataset), so the
# sweep can generate independent simulation replicates.
# =============================================================================

# NB: MASS (mvrnorm, rnegbin) is used via MASS:: prefixes rather than
# library(MASS) so it does not mask dplyr::select downstream.
suppressMessages(library(tidyverse))

# -----------------------------------------------------------------------------
# Default vS1 GP truth (shared lx/ly across species; differing marginal SD).
# Exposed as an argument so a future sweep can stress shorter length scales.
# -----------------------------------------------------------------------------
vS1_default_gp_params <- function() {
  sp_names <- c("Merluccius_productus",
                "Megaptera_novaeangliae",
                "Lagenorhynchus_obliquidens")
  list(
    hake     = list(name = sp_names[1], sigma = 2.0, lx = 50, ly = 300,
                    mu = log(15),   zsample_pref = 0.0),
    humpback = list(name = sp_names[2], sigma = 1.0, lx = 50, ly = 300,
                    mu = log(0.01), zsample_pref = 0.0),
    pwsd     = list(name = sp_names[3], sigma = 1.3, lx = 50, ly = 300,
                    mu = log(0.02), zsample_pref = 0.0)
  )
}

# -----------------------------------------------------------------------------
# simulate_whale_edna_vS1(): one vS1 simulation replicate.
# Returns the same nested list structure that 01 saves to disk, so
# format_stan_data_vS1() can consume it directly.
# -----------------------------------------------------------------------------
simulate_whale_edna_vS1 <- function(seed = 1L,
                                    n_stations = 250L,
                                    gp_params  = vS1_default_gp_params()) {
  set.seed(seed)

  n_species <- length(gp_params)
  sp_names  <- vapply(gp_params, `[[`, character(1), "name")
  sp_common <- c("Pacific hake", "Humpback whale", "Pacific white-sided dolphin")
  has_qpcr  <- c(TRUE, FALSE, FALSE)
  n_qpcr_rep <- 3L; n_mb_rep <- 3L
  prop_reduced_qpcr_rep <- 0; prop_reduced_mb_rep <- 0

  X_min <- 100000; X_max <- 600000
  Y_min <- 4180000; Y_max <- 5450000
  X_km_max <- (X_max - X_min) / 1000
  Y_km_max <- (Y_max - Y_min) / 1000

  conv_factor    <- c(hake = 10, humpback = 200, pwsd = 110, junk = 10)
  sample_depths  <- 0; n_sample_depth <- 1L
  vol_filtered   <- 2.5; vol_aliquot <- 2

  mb_reads_p_tight       <- 0.80
  mb_reads_tight_meanlog <- log(75000); mb_reads_tight_sdlog <- 0.25
  mb_reads_wide_meanlog  <- log(40000); mb_reads_wide_sdlog  <- 1.00
  mb_reads_min <- 1000L; mb_reads_max <- 250000L

  # --- stations + samples (surface only) ---
  stations <- tibble(
    station = seq_len(n_stations),
    X = runif(n_stations, 0, X_km_max),
    Y = runif(n_stations, 0, Y_km_max),
    X_utm = X_min + X * 1000,
    Y_utm = Y_min + Y * 1000
  )
  samples <- expand.grid(station = 1:n_stations, depth_idx = 1:n_sample_depth) %>%
    left_join(stations, by = "station") %>%
    mutate(Z_sample = sample_depths[depth_idx], sample_id = row_number())
  N <- nrow(samples)
  coords_gp <- as.matrix(samples[, c("X", "Y")])

  # --- 2-D anisotropic zero-mean GP per species ---
  aniso_cov <- function(coords, sigma, lx, ly) {
    n <- nrow(coords); K <- matrix(0, n, n)
    for (i in seq_len(n)) for (j in i:n) {
      d2 <- ((coords[i,1]-coords[j,1])/lx)^2 + ((coords[i,2]-coords[j,2])/ly)^2
      K[i,j] <- K[j,i] <- sigma^2 * exp(-0.5 * d2)
    }
    K + diag(1e-6, n)
  }
  lambda_true_si <- matrix(NA, N, n_species)
  gp_field_si    <- matrix(NA, N, n_species)
  for (s in seq_len(n_species)) {
    p  <- gp_params[[s]]
    K  <- aniso_cov(coords_gp, p$sigma, p$lx, p$ly)
    fs <- as.vector(MASS::mvrnorm(1, mu = rep(0.0, N), Sigma = K))
    gp_field_si[, s]    <- fs
    lambda_true_si[, s] <- exp(p$mu + fs)
  }

  # --- surface-only water-column effect (exp(0) = 1) ---
  zsample_effect <- matrix(1.0, N, n_species)

  # --- bottle copies ---
  C_obs_si <- matrix(NA_integer_, N, n_species)
  for (s in seq_len(n_species)) {
    C_obs_si[, s] <- MASS::rnegbin(
      N, mu = conv_factor[s] * lambda_true_si[, s] * zsample_effect[, s] * vol_filtered,
      theta = 10)
  }

  # --- replicate counts (full design here) ---
  reduced_qpcr <- runif(N) < prop_reduced_qpcr_rep
  n_qpcr_rep_i <- ifelse(reduced_qpcr, sample.int(2L, N, TRUE), as.integer(n_qpcr_rep))
  reduced_mb   <- runif(N) < prop_reduced_mb_rep
  n_mb_rep_i   <- ifelse(reduced_mb, sample.int(2L, N, TRUE), as.integer(n_mb_rep))
  N_qpcr_long  <- sum(n_qpcr_rep_i); N_mb_long <- sum(n_mb_rep_i)
  qpcr_sample_idx <- rep(seq_len(N), times = n_qpcr_rep_i)
  mb_sample_idx   <- rep(seq_len(N), times = n_mb_rep_i)

  # --- aliquot copies ---
  qpcr_copies     <- rbinom(N_qpcr_long, size = C_obs_si[qpcr_sample_idx, 1], prob = vol_aliquot/100)
  qpcr_exp_copies <- C_obs_si[qpcr_sample_idx, 1] * (vol_aliquot/100)
  mb_copies <- matrix(NA_integer_, N_mb_long, n_species)
  for (s in seq_len(n_species))
    mb_copies[, s] <- rbinom(N_mb_long, size = C_obs_si[mb_sample_idx, s], prob = vol_aliquot/100)

  # --- qPCR hurdle (hake only) ---
  qpcr_p <- list(kappa = 0.85, alpha_ct = 38.0, beta_ct = 1.44,
                 gamma_0 = 0.30, gamma_1 = -0.30, sigma_0 = 0.23)
  p_det_hake  <- 1 - exp(-qpcr_p$kappa * qpcr_exp_copies)
  qpcr_detect <- rbinom(N_qpcr_long, 1, p_det_hake)
  qpcr_p$sigma_ct <- sqrt(qpcr_p$sigma_0^2 +
                          exp(2*(qpcr_p$gamma_0 + qpcr_p$gamma_1 * log(qpcr_exp_copies))))
  qpcr_ct <- ifelse(qpcr_detect == 1L,
                    qpcr_p$alpha_ct - qpcr_p$beta_ct * log(pmax(qpcr_exp_copies,1)) +
                      rnorm(N_qpcr_long, 0, qpcr_p$sigma_ct),
                    NA_real_)

  # --- metabarcoding: junk + target multinomial ---
  junk_mean <- 50; junk_theta <- 1000
  junk_exp_copies <- MASS::rnegbin(nrow(samples), mu = junk_mean, theta = junk_theta)
  C_obs_si_junk   <- junk_exp_copies * 50
  # junk drawn once per sample, then expanded to aliquot long form so it is
  # CONSTANT within a sample and aligned with the per-aliquot target draws.
  mb_copies_junk      <- rbinom(N, size = C_obs_si_junk, prob = vol_aliquot/100)
  mb_copies_junk_long <- mb_copies_junk[mb_sample_idx]

  component_tight <- rbinom(N_mb_long, 1, mb_reads_p_tight)
  read_depth <- ifelse(component_tight == 1L,
                       rlnorm(N_mb_long, mb_reads_tight_meanlog, mb_reads_tight_sdlog),
                       rlnorm(N_mb_long, mb_reads_wide_meanlog,  mb_reads_wide_sdlog))
  read_depth <- pmin(pmax(as.integer(round(read_depth)), mb_reads_min), mb_reads_max)

  sample_total <- rowSums(mb_copies) + mb_copies_junk_long
  pi_edna <- mb_copies           / sample_total
  pi_junk <- mb_copies_junk_long / sample_total
  pi_all  <- cbind(pi_edna, pi_junk)

  mb_reads <- t(vapply(seq_len(N_mb_long),
                       function(i) rmultinom(1, size = read_depth[i], prob = pi_all[i, ])[, 1],
                       integer(n_species + 1)))
  storage.mode(mb_reads) <- "integer"
  mb_total <- as.integer(rowSums(mb_reads))

  list(
    meta = list(n_species = n_species, sp_names = sp_names, sp_common = sp_common,
                has_qpcr = has_qpcr, conv_factor = conv_factor,
                vol_filtered = vol_filtered, vol_aliquot = vol_aliquot,
                N = N, N_qpcr_long = N_qpcr_long, N_mb_long = N_mb_long,
                X_km_max = X_km_max, Y_km_max = Y_km_max, seed = seed),
    design = list(stations = stations, samples = samples,
                  n_qpcr_rep_i = n_qpcr_rep_i, n_mb_rep_i = n_mb_rep_i),
    truth = list(gp_field_si = gp_field_si, lambda_true_si = lambda_true_si,
                 C_obs_si = C_obs_si, zsample_effect = zsample_effect,
                 gp_params = gp_params, qpcr_params = qpcr_p),
    observed = list(qpcr_sample_idx = qpcr_sample_idx, qpcr_detect = qpcr_detect,
                    qpcr_ct = qpcr_ct, mb_sample_idx = mb_sample_idx,
                    mb_reads = mb_reads, mb_total = mb_total)
  )
}

# -----------------------------------------------------------------------------
# format_stan_data_vS1(): build the Stan data list for a given basis count.
# Mirrors 03. prediction_n keeps a (tiny by default) prediction grid so the
# model's generated-quantities block stays valid; the sweep does not use the
# predictions, so a small grid keeps fits fast and CSVs small.
# -----------------------------------------------------------------------------
format_stan_data_vS1 <- function(sim,
                                 HSGP_M = c(16L, 8L),
                                 HSGP_C = c(1.5, 1.5),
                                 prediction_n = 4L) {
  samples        <- sim$design$samples
  zsample_effect <- sim$truth$zsample_effect
  N            <- sim$meta$N
  S            <- sim$meta$n_species
  vol_aliquot  <- sim$meta$vol_aliquot
  vol_filtered <- sim$meta$vol_filtered
  conv_factor  <- sim$meta$conv_factor
  X_km_max     <- sim$meta$X_km_max
  Y_km_max     <- sim$meta$Y_km_max

  # water-column log-offset
  log_zsample_effect <- log(as.matrix(zsample_effect))
  log_zsample_effect[!is.finite(log_zsample_effect)] <- -10.0

  # coordinates normalised to [-1, 1]
  coords_raw   <- cbind(as.numeric(samples[["X"]]), as.numeric(samples[["Y"]]))
  coord_centre <- c(X_km_max/2, Y_km_max/2)
  coord_scale  <- c(X_km_max/2, Y_km_max/2)
  coords_norm  <- sweep(sweep(coords_raw, 2, coord_centre, "-"), 2, coord_scale, "/")
  stopifnot(max(abs(coords_norm)) < 1.5)

  # small prediction grid (unused downstream but keeps GQ valid)
  npc <- max(2L, round(sqrt(prediction_n)))
  pred_points <- expand.grid(x = seq(0, X_km_max, length.out = npc),
                             y = seq(0, Y_km_max, length.out = npc))
  N_pred <- nrow(pred_points)
  pred_coords_norm <- sweep(sweep(as.matrix(pred_points), 2, coord_centre, "-"),
                            2, coord_scale, "/")

  # qPCR passthrough
  qpcr_sample_idx <- as.integer(sim$observed$qpcr_sample_idx)
  qpcr_detect_vec <- as.integer(sim$observed$qpcr_detect)
  qpcr_ct_vec     <- as.numeric(sim$observed$qpcr_ct); qpcr_ct_vec[is.na(qpcr_ct_vec)] <- 0.0

  # MB passthrough, drop zero-read rows
  mb_reads_full <- sim$observed$mb_reads; storage.mode(mb_reads_full) <- "integer"
  mb_total_full <- as.integer(sim$observed$mb_total)
  keep <- mb_total_full > 0
  mb_reads_long <- mb_reads_full[keep, , drop = FALSE]
  mb_sample_idx <- as.integer(sim$observed$mb_sample_idx)[keep]
  mb_total_vec  <- mb_total_full[keep]
  N_mb_long     <- nrow(mb_reads_long)

  # basis-function index grid
  D1 <- length(HSGP_M)
  INDICES <- as.matrix(do.call(expand_grid, lapply(HSGP_M, seq_len)))

  stan_data <- list(
    N = as.integer(N), N_pred = as.integer(N_pred), S = as.integer(S),
    M = as.integer(prod(HSGP_M)), D1 = as.integer(D1), INDICES = INDICES,
    coords = coords_norm, coord_scale = coord_scale, pred_coords = pred_coords_norm,
    L_hsgp = HSGP_C, m_hsgp = as.integer(HSGP_M),
    log_zsample_effect = log_zsample_effect,
    log_conv_factor = as.numeric(log(conv_factor)),
    vol_aliquot = vol_aliquot, log_vol_filtered = log(vol_filtered),
    N_qpcr_long = length(qpcr_sample_idx), qpcr_sample_idx = qpcr_sample_idx,
    qpcr_detect = qpcr_detect_vec, qpcr_ct = qpcr_ct_vec,
    N_mb_long = N_mb_long, mb_sample_idx = mb_sample_idx,
    mb_reads = mb_reads_long, mb_total = mb_total_vec,
    alpha_ct = sim$truth$qpcr_params$alpha_ct, beta_ct = sim$truth$qpcr_params$beta_ct,
    kappa = sim$truth$qpcr_params$kappa, gamma0_ct = sim$truth$qpcr_params$gamma_0,
    gamma1_ct = sim$truth$qpcr_params$gamma_1, sigma0_ct = sim$truth$qpcr_params$sigma_0,
    prior_mu_sp_mu = 2.0, prior_mu_sp_sig = 1.5,
    prior_gp_sigma_shape = 8.0, prior_gp_sigma_rate = 4.0,
    prior_gp_raw_alpha = 10, prior_gp_raw_beta = 16,
    prior_beta0_phi_mu = 0.0, prior_beta0_phi_sig = 1.0,
    prior_gamma0_phi_mu = 2.0, prior_gamma0_phi_sig = 1.0,
    prior_gamma1_phi_mu = 1.0, prior_gamma1_phi_sig = 0.5
  )
  list(stan_data = stan_data, pred_points = pred_points)
}
