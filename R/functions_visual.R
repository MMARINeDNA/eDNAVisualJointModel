# =============================================================================
# functions_visual.R
#
# Simulation + Stan-data formatting for the spatial line-transect (visual)
# model, stan/hsgp_visual.stan. Run from the repo root.
#
#   sim <- simulate_visual(seed = 1)                    # one truth, both species
#   sd  <- format_stan_data_visual(sim, "humpback")     # Stan data for one species
#
# The latent field is HSGP4eDNA's `bathysp` truth, drawn by its
# simulate_field_bathysp():
#   log lambda_s(x) = mu_s + f_s(X, Y) + B(Z_bathy) . beta_s     (animals / km^2)
# (default_gp_params(), bathy_spline_basis(), BATHY_BETA_TRUE), so visual and
# eDNA data can be generated from ONE field (R/functions_joint.R).
#
# The detection / group-size simulation is distance/00_distance_v4.1.R's:
# half-normal detection (PWSD with a group-size covariate), group sizes
# resampled from the empirical pools in data/grpsz/. With sea_state = TRUE each
# segment also gets a Beaufort sea state with a spatial trend, and detection
# scales as exp(beta_bft * (Beaufort - BFT_CENTRE)). With a finite
# `nb_size_truth` per species, groups per segment are negative binomial
# (clustered sightings); the default (Inf) is Poisson.
# =============================================================================

HSGP4EDNA_DIR <- "external/HSGP4eDNA"
source(file.path(HSGP4EDNA_DIR, "R", "functions.R"))   # default_gp_params(),
                                                       # hsgp_basis_rule(),
                                                       # bathy_spline_basis(),
                                                       # gp_prior_from_design()

# Study domain (km from the south-west corner, UTM 10N): HSGP4eDNA's
# STUDY_DOMAIN, the MURI study area north of 35°N (690 x 1770 km). Shared by
# the eDNA, visual and joint simulators and the real-data loader.
LT_DOMAIN <- list(X_km_max = STUDY_DOMAIN$X_km_max, Y_km_max = STUDY_DOMAIN$Y_km_max)

# Bottom depth (m), drawn independently of (X, Y) as in simulate_bathysp(), so
# the spline effect is identifiable separately from the spatial GP.
draw_bathy <- function(n) pmax(50, 50 + 3150 * rbeta(n, 2, 2))

# Sea state. Detection scales are defined at Beaufort BFT_CENTRE, so sigma_det
# and the model's sigma_c are the scale in Beaufort 3.
BFT_CENTRE <- 3

# Segment Beaufort with a spatial trend: rougher to the north (large Y) and
# offshore (small X, the west of the domain), plus segment noise; on [0, 6].
# Mean ~3, trend +/- 1.75. Ignoring it in the model then biases the field
# along the trend, which is what the sea-state validation checks.
draw_sea_state <- function(segments, domain = LT_DOMAIN, sd = 0.5) {
  b <- BFT_CENTRE + 1.0 * (2 * segments$Y / domain$Y_km_max - 1) +
    0.75 * (1 - 2 * segments$X / domain$X_km_max) + rnorm(nrow(segments), 0, sd)
  pmin(pmax(b, 0), 6)
}

# -----------------------------------------------------------------------------
# species_field(): HSGP4eDNA's simulate_field_bathysp() for the named species,
# with species-named columns. One GP draw per species, in `species` order.
# -----------------------------------------------------------------------------
species_field <- function(locs, species = c("humpback", "pwsd"),
                          gp_params = default_gp_params(), df = 4L,
                          spline_setup = NULL, bathy_scale = 1.0) {
  f <- simulate_field_bathysp(locs, gp_params[species], df = df,
                              spline_setup = spline_setup, bathy_scale = bathy_scale)
  nm <- function(m) { colnames(m) <- species; m }
  list(log_lambda = nm(f$log_lambda_si), gp_field = nm(f$gp_field_si),
       bathy_effect = nm(f$bathy_effect_si), spline_setup = f$bathy_spline,
       gp_params = gp_params[species],
       beta_bathy_true = `rownames<-`(f$beta_bathy_true, species))
}

# -----------------------------------------------------------------------------
# lt_design(): systematic E-W transects split into segments (as v4.1).
# -----------------------------------------------------------------------------
lt_design <- function(n_transects = 35L, seg_length = 10, domain = LT_DOMAIN) {
  n_seg_per <- floor(domain$X_km_max / seg_length)
  transect_Y <- seq(domain$Y_km_max / (2 * n_transects),
                    domain$Y_km_max - domain$Y_km_max / (2 * n_transects),
                    length.out = n_transects)
  seg <- expand.grid(transect_id = seq_len(n_transects), seg_idx = seq_len(n_seg_per))
  data.frame(seg_id = seq_len(nrow(seg)), transect_id = seg$transect_id,
             X = (seg$seg_idx - 0.5) * seg_length, Y = transect_Y[seg$transect_id],
             seg_l = seg_length)
}

# -----------------------------------------------------------------------------
# lt_species_params(): detection, group-size and prior settings per species,
# from distance/00_distance_v4.1.R. Density/GP truths are NOT here - they come
# from default_gp_params() so they are shared with the eDNA model.
# -----------------------------------------------------------------------------
lt_species_params <- function(grpsz_dir = "data/grpsz") {
  pool_h <- as.integer(readRDS(file.path(grpsz_dir, "humpback.rds")))
  pool_p <- as.integer(readRDS(file.path(grpsz_dir, "pwsd.rds")))
  mean_h <- mean(pool_h); mean_p <- mean(pool_p)
  cv2_p  <- (sd(pool_p) / mean_p)^2
  pwsd_log_mu0    <- log(mean_p) - 0.5 * log(1 + cv2_p)   # log-normal MoM
  pwsd_log_sigma0 <- sqrt(log(1 + cv2_p))
  list(
    humpback = list(
      common_name = "Humpback whale", group_size_pool = pool_h, mean_group_size = mean_h,
      sigma_det = 2.5, use_size_covar = 0L, beta_size_truth = 0.0, s_centre = mean_h,
      beta_bft_truth = -0.08, nb_size_truth = Inf,
      det_key = 0L, hr_b_truth = 3, log_b_prior_mean = log(3), log_b_prior_sd = 0.7,
      model_group_dist = 0L, S_max = 50L,
      log_sigma_prior_mean = log(2.5), log_sigma_prior_sd = 0.6,
      beta_size_prior_mean = 0.0, beta_size_prior_sd = 0.05,
      beta_bft_prior_mean = -0.1, beta_bft_prior_sd = 0.2,
      mu_s_prior_shape = 4, mu_s_prior_rate = 4 / mean_h,
      phi_s_prior_shape = 1, phi_s_prior_rate = 0.1,
      mu_log_prior_mean = log(mean_h), mu_log_prior_sd = 1.0,
      sigma_log_prior_shape = 2, sigma_log_prior_rate = 2
    ),
    pwsd = list(
      common_name = "Pacific white-sided dolphin", group_size_pool = pool_p, mean_group_size = mean_p,
      sigma_det = 1.5, use_size_covar = 1L, beta_size_truth = 0.01, s_centre = mean_p,
      beta_bft_truth = -0.15, nb_size_truth = Inf,
      det_key = 0L, hr_b_truth = 3, log_b_prior_mean = log(3), log_b_prior_sd = 0.7,
      model_group_dist = 1L, S_max = 1000L,
      log_sigma_prior_mean = log(1.5), log_sigma_prior_sd = 0.6,
      beta_size_prior_mean = 0.0, beta_size_prior_sd = 0.05,
      beta_bft_prior_mean = -0.1, beta_bft_prior_sd = 0.2,
      mu_s_prior_shape = 4, mu_s_prior_rate = 4 / mean_p,
      phi_s_prior_shape = 1, phi_s_prior_rate = 0.1,
      mu_log_prior_mean = pwsd_log_mu0, mu_log_prior_sd = 0.5,
      sigma_log_prior_shape = 4, sigma_log_prior_rate = 4 / pwsd_log_sigma0
    )
  )
}

# -----------------------------------------------------------------------------
# simulate_lt_sightings(): groups in each segment's strip ~ Poisson, uniform
# perpendicular distance, half-normal detection (v4.1 logic). If segments has
# a `bft` column, each segment's detection scale is multiplied by
# exp(beta_bft_truth * (bft - BFT_CENTRE)); without it the draws are v4.1's.
# Clustering: with a finite p$nb_size_truth (k), each segment's expected
# groups are multiplied by a gamma(k, k) draw (mean 1, variance 1/k), so
# groups per segment are negative binomial with size k; with k = Inf (or
# unset) no gamma draws are made and the Poisson draws are unchanged.
# lambda_animals: animals / km^2 per segment.
# -----------------------------------------------------------------------------
simulate_lt_sightings <- function(segments, lambda_animals, p, w) {
  lambda_groups <- lambda_animals / p$mean_group_size
  mean_strip <- lambda_groups * 2 * w * segments$seg_l
  k <- if (is.null(p$nb_size_truth)) Inf else p$nb_size_truth
  if (is.finite(k)) mean_strip <- mean_strip * rgamma(nrow(segments), shape = k, rate = k)
  n_strip <- rpois(nrow(segments), mean_strip)
  rows <- vector("list", nrow(segments))
  for (j in which(n_strip > 0)) {
    x_true  <- runif(n_strip[j], 0, w)
    s_strip <- sample(p$group_size_pool, n_strip[j], replace = TRUE)
    sigma_g <- if (p$use_size_covar == 1L)
      p$sigma_det * exp(p$beta_size_truth * (s_strip - p$s_centre)) else
      rep(p$sigma_det, n_strip[j])
    if (!is.null(segments$bft))
      sigma_g <- sigma_g * exp(p$beta_bft_truth * (segments$bft[j] - BFT_CENTRE))
    g <- if (isTRUE(p$det_key == 1L))
      -expm1(-(sigma_g / pmax(x_true, 1e-6))^p$hr_b_truth) else
      exp(-x_true^2 / (2 * sigma_g^2))
    keep <- runif(n_strip[j]) < g
    if (any(keep))
      rows[[j]] <- data.frame(seg_id = segments$seg_id[j], distance = x_true[keep],
                              size = as.integer(s_strip[keep]))
  }
  obs <- do.call(rbind, rows)
  if (is.null(obs)) obs <- data.frame(seg_id = integer(), distance = numeric(), size = integer())
  list(obs = obs, seg_count = tabulate(obs$seg_id, nbins = nrow(segments)),
       n_groups_in_strip = n_strip)
}

# -----------------------------------------------------------------------------
# simulate_visual(): one field over the segment design + sightings per species.
# sea_state = TRUE adds segment Beaufort (draw_sea_state(), drawn after the
# field, so the field matches the sea_state = FALSE draw for the same seed).
# segments: a design to use instead of lt_design() + draw_bathy(), e.g. the
#   real effort (data/processed/muri_visual.rds design$segments); it needs
#   seg_id (1..n), transect_id, X, Y (km), seg_l (km) and Z_bathy (m).
# nb_size: named NB size per species (clustered sightings), overriding
#   sp_params' nb_size_truth; Inf = Poisson.
# -----------------------------------------------------------------------------
simulate_visual <- function(seed = 1L, species = c("humpback", "pwsd"),
                            n_transects = 35L, seg_length = 10, w = 5.0, df = 4L,
                            gp_params = default_gp_params(),
                            sp_params = lt_species_params(),
                            sea_state = FALSE, segments = NULL, nb_size = NULL) {
  set.seed(seed)
  if (is.null(segments)) {
    segments <- lt_design(n_transects, seg_length)
    segments$Z_bathy <- draw_bathy(nrow(segments))
  } else {
    segments <- as.data.frame(segments)[, c("seg_id", "transect_id", "X", "Y", "seg_l", "Z_bathy")]
    stopifnot(identical(segments$seg_id, seq_len(nrow(segments))), !anyNA(segments))
  }
  for (sp in intersect(names(nb_size), species)) sp_params[[sp]]$nb_size_truth <- nb_size[[sp]]
  field <- species_field(segments, species, gp_params = gp_params, df = df)
  if (sea_state) segments$bft <- draw_sea_state(segments)
  observed <- lapply(setNames(species, species), function(sp)
    simulate_lt_sightings(segments, exp(field$log_lambda[, sp]), sp_params[[sp]], w))
  list(meta = list(seed = seed, species = species, w = w, domain = LT_DOMAIN),
       design = list(segments = segments),
       truth = c(field, list(sp_params = sp_params[species])),
       observed = observed)
}

# -----------------------------------------------------------------------------
# format_stan_data_visual(): Stan data for one species.
#   HSGP_M    basis per axis; NULL = sized from the design-based prior
#             (gp_prior_from_design(), within the M_max basis budget)
#   use_gp    0 = no spatial field
#   use_bathy FALSE = no bottom-depth spline
#   use_bft   TRUE = sea-state covariate on detection (default: the species'
#             `use_bft` setting in sp_params if present, else whenever the
#             segments carry a `bft` column). Beaufort is binned to bft_bin
#             steps; the size-averaged ESW is computed once per level.
#   use_nb    1 = negative-binomial segment counts (default: the species'
#             `use_nb` in sp_params, else 0 = Poisson); nb_inv_sqrt_sd is the
#             half-normal prior sd on 1 / sqrt(phi_nb)
#   Detection key: sp_params$det_key (0 half-normal, 1 hazard-rate).
#   Truncation: sim$meta$w, one value or named per species. Detections beyond
#   it are dropped here and seg_count is rebuilt, so objects may carry all
#   sightings while each species gets its own w.
# Coordinates are normalised by the DOMAIN extents (not the data range), so the
# same normalisation can serve eDNA stations and LT segments in the joint model.
# -----------------------------------------------------------------------------
format_stan_data_visual <- function(sim, sp, HSGP_M = NULL, HSGP_C = c(1.5, 1.5),
                                    use_gp = 1L, use_bathy = TRUE,
                                    prior_mu_sp = c(-5, 2),
                                    M_max = GP_BASIS_BUDGET,
                                    use_bft = {
                                      ub <- sim$truth$sp_params[[sp]]$use_bft
                                      if (is.null(ub)) !is.null(sim$design$segments$bft) else isTRUE(ub)
                                    },
                                    bft_bin = 0.5, bft_centre = BFT_CENTRE,
                                    use_nb = as.integer(isTRUE(sim$truth$sp_params[[sp]]$use_nb == 1L)),
                                    nb_inv_sqrt_sd = 1.0) {
  seg <- sim$design$segments
  p   <- sim$truth$sp_params[[sp]]
  dom <- sim$meta$domain
  `%||%` <- function(a, b) if (is.null(a)) b else a
  w_all <- sim$meta$w
  w <- if (!is.null(names(w_all)) && sp %in% names(w_all)) w_all[[sp]] else as.numeric(w_all)[1]
  ob_all <- sim$observed[[sp]]
  o <- ob_all$obs[ob_all$obs$distance <= w, , drop = FALSE]
  ob <- list(obs = o, seg_count = tabulate(o$seg_id, nbins = nrow(seg)))
  coord_centre <- c(dom$X_km_max, dom$Y_km_max) / 2
  coord_scale  <- coord_centre
  coords <- sweep(sweep(as.matrix(seg[, c("X", "Y")]), 2, coord_centre, "-"), 2, coord_scale, "/")
  # Design-based GP priors (HSGP4eDNA gp_prior_from_design()): segment locations only - no
  # truth and no detections are used. The basis is
  # sized to the prior's lower length-scale bound unless HSGP_M is given.
  gp_prior <- gp_prior_from_design(seg[, c("X", "Y")],
                                   domain_range = c(dom$X_km_max, dom$Y_km_max),
                                   coord_scale = coord_scale,
                                   M_max = M_max, c = HSGP_C[1])
  if (is.null(HSGP_M)) HSGP_M <- gp_prior$HSGP_M
  INDICES <- as.matrix(do.call(tidyr::expand_grid, lapply(HSGP_M, seq_len)))
  B <- if (use_bathy) bathy_spline_basis(seg$Z_bathy, setup = sim$truth$spline_setup)$B else
    matrix(0, nrow(seg), 0)
  # Sea-state levels: one level (and no covariate) unless use_bft
  if (use_bft) {
    stopifnot(!is.null(seg$bft), !anyNA(seg$bft))
    b_bin    <- round(seg$bft / bft_bin) * bft_bin
    levels   <- sort(unique(b_bin))
    seg_bft  <- match(b_bin, levels)
  } else {
    levels <- bft_centre; seg_bft <- rep(1L, nrow(seg))
  }
  det_bft <- seg_bft[match(ob$obs$seg_id, seg$seg_id)]

  sd <- list(
    n = nrow(ob$obs), x = ob$obs$distance, s = as.integer(ob$obs$size), w = w,
    det_key = as.integer(p$det_key %||% 0L),
    log_b_prior_mean = p$log_b_prior_mean %||% log(3), log_b_prior_sd = p$log_b_prior_sd %||% 0.7,
    log_sigma_prior_mean = p$log_sigma_prior_mean, log_sigma_prior_sd = p$log_sigma_prior_sd,
    use_size_covar = p$use_size_covar, s_centre = p$s_centre,
    beta_size_prior_mean = p$beta_size_prior_mean, beta_size_prior_sd = p$beta_size_prior_sd,
    use_bft_covar = as.integer(use_bft), n_bft = length(levels),
    bft_level = as.array(levels - bft_centre),
    beta_bft_prior_mean = if (is.null(p$beta_bft_prior_mean)) -0.1 else p$beta_bft_prior_mean,
    beta_bft_prior_sd   = if (is.null(p$beta_bft_prior_sd))   0.2 else p$beta_bft_prior_sd,
    seg_bft = as.integer(seg_bft), det_bft = as.array(as.integer(det_bft)),
    use_nb = as.integer(use_nb), prior_nb_inv_sqrt_sd = nb_inv_sqrt_sd,
    S_max = p$S_max, group_size_dist = p$model_group_dist,
    mu_s_prior_shape = p$mu_s_prior_shape, mu_s_prior_rate = p$mu_s_prior_rate,
    phi_s_prior_shape = p$phi_s_prior_shape, phi_s_prior_rate = p$phi_s_prior_rate,
    mu_log_prior_mean = p$mu_log_prior_mean, mu_log_prior_sd = p$mu_log_prior_sd,
    sigma_log_prior_shape = p$sigma_log_prior_shape, sigma_log_prior_rate = p$sigma_log_prior_rate,
    n_seg = nrow(seg), seg_count = as.integer(ob$seg_count), seg_l = as.numeric(seg$seg_l),
    use_gp = as.integer(use_gp), D1 = 2L, M = as.integer(prod(HSGP_M)), INDICES = INDICES,
    coords = coords, coord_scale = coord_scale, L_hsgp = HSGP_C,
    K_bathy = ncol(B), B_bathy = B,
    N_pred = 0L, pred_coords = matrix(0, 0, 2), B_bathy_pred = matrix(0, 0, ncol(B)),
    # Field priors: gp_sigma and gp_l from the design rule; mu_sp weakly
    # informative on log animals / km^2 (-5 +/- 2 -> 1e-4 .. 0.4);
    # beta_bathy as HSGP4eDNA.
    prior_mu_sp_mu = prior_mu_sp[1], prior_mu_sp_sig = prior_mu_sp[2],
    prior_gp_sigma_shape = gp_prior$prior_gp_sigma_shape,
    prior_gp_sigma_rate  = gp_prior$prior_gp_sigma_rate,
    prior_gp_l_shape = gp_prior$prior_gp_l_shape,
    prior_gp_l_scale = gp_prior$prior_gp_l_scale,
    prior_beta_bathy_sig = 2.0
  )
  attr(sd, "gp_prior") <- gp_prior   # provenance (an attribute, so cmdstanr ignores it)
  sd
}
