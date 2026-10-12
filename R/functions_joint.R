# =============================================================================
# functions_joint.R
#
# Joint simulator: ONE latent animal-density field per species, observed by
# both eDNA sampling (stations) and visual line-transect surveys (segments).
#
#   sim <- simulate_joint(seed = 1)
#   sim$edna     # shaped like HSGP4eDNA's simulate_bathysp() output
#   sim$visual   # shaped like simulate_visual() output
#
# The field is HSGP4eDNA's `bathysp` truth (simulate_field_bathysp()), drawn
# jointly over the union of station and segment locations so the two data
# sources see the same realisation:
#   log lambda_s(x) = mu_s + f_s(X, Y) + B(Z_bathy) . beta_s     (animals / km^2)
# The bottom-depth spline basis is built from all locations (stations and
# segments) so both sources share one h_s(Z_bathy).
#
# Observation models are unchanged from the single-source simulators:
#   eDNA   simulate_edna_obs() (HSGP4eDNA): qPCR on hake, junk-background
#          metabarcoding on all species
#   visual simulate_lt_sightings(): half-normal detection, empirical group
#          sizes, for the cetaceans (lt_species)
# Because each sub-object has its single-source shape, every eDNA-only and
# visual-only fitting tool works on it unchanged.
# =============================================================================

source("R/functions_visual.R")

simulate_joint <- function(seed = 1L, n_stations = 200L,
                           sample_depths = 0, zsample_pref = NULL,
                           n_transects = 35L, seg_length = 10, w = 5.0, df = 4L,
                           gp_params  = default_gp_params(),
                           lt_species = c("humpback", "pwsd"),
                           sp_params  = lt_species_params()) {
  set.seed(seed)
  dom <- LT_DOMAIN
  S <- length(gp_params)
  sp_common <- c("Pacific hake", "Humpback whale", "Pacific white-sided dolphin")
  conv_factor <- c(hake = 10, humpback = 200, pwsd = 110, junk = 10)
  vol_filtered <- 2.5; vol_aliquot <- 2

  # ---- Designs -----------------------------------------------------------
  # eDNA stations as in simulate_bathysp(); samples = station x sample depth,
  # dropping samples below the sea floor.
  stations <- data.frame(station = seq_len(n_stations),
                         X = runif(n_stations, 0, dom$X_km_max),
                         Y = runif(n_stations, 0, dom$Y_km_max))
  stations$Z_bathy <- draw_bathy(n_stations)
  samples <- tidyr::expand_grid(stations, Z_sample = sample_depths)
  samples <- samples[samples$Z_sample <= samples$Z_bathy, ]
  samples$depth_idx <- match(samples$Z_sample, sample_depths)
  samples$sample_id <- seq_len(nrow(samples))

  segments <- lt_design(n_transects, seg_length, dom)
  segments$Z_bathy <- draw_bathy(nrow(segments))

  # ---- One field over all locations (stations first, then segments) --------
  locs <- rbind(stations[, c("X", "Y", "Z_bathy")], segments[, c("X", "Y", "Z_bathy")])
  locs$source <- rep(c("station", "segment"), c(n_stations, nrow(segments)))
  field <- simulate_field_bathysp(locs, gp_params, df = df)
  st  <- seq_len(n_stations)
  sgi <- n_stations + seq_len(nrow(segments))

  # ---- eDNA observations ---------------------------------------------------
  rows <- st[samples$station]                       # sample -> field row
  lambda_si <- exp(field$log_lambda_si[rows, , drop = FALSE])
  zse <- zsample_effect_matrix(samples$Z_sample, S, sample_depths, zsample_pref)
  eo  <- simulate_edna_obs(lambda_si, zse, conv_factor = conv_factor,
                           vol_filtered = vol_filtered, vol_aliquot = vol_aliquot)
  edna <- list(
    meta = list(n_species = S, sp_names = vapply(gp_params, `[[`, character(1), "name"),
                sp_common = sp_common, conv_factor = conv_factor,
                vol_filtered = vol_filtered, vol_aliquot = vol_aliquot,
                N = nrow(samples), N_qpcr_long = eo$N_qpcr_long, N_mb_long = eo$N_mb_long,
                X_km_max = dom$X_km_max, Y_km_max = dom$Y_km_max, seed = seed),
    design = list(stations = stations, samples = samples,
                  n_qpcr_rep_i = eo$n_qpcr_rep_i, n_mb_rep_i = eo$n_mb_rep_i),
    truth = list(gp_field_si = field$gp_field_si[rows, , drop = FALSE],
                 bathy_effect_si = field$bathy_effect_si[rows, , drop = FALSE],
                 lambda_true_si = lambda_si, C_obs_si = eo$C_obs_si,
                 zsample_effect = zse, gp_params = gp_params, qpcr_params = eo$qpcr_params,
                 beta_bathy_true = field$beta_bathy_true, bathy_spline = field$bathy_spline),
    observed = eo$observed)

  # ---- Visual observations -------------------------------------------------
  pick <- function(m) { m <- m[sgi, match(lt_species, names(gp_params)), drop = FALSE]
                        colnames(m) <- lt_species; m }
  log_lambda_seg <- pick(field$log_lambda_si)
  observed_lt <- lapply(setNames(lt_species, lt_species), function(sp)
    simulate_lt_sightings(segments, exp(log_lambda_seg[, sp]), sp_params[[sp]], w))
  visual <- list(
    meta = list(seed = seed, species = lt_species, w = w, domain = dom),
    design = list(segments = segments),
    truth = list(log_lambda = log_lambda_seg, gp_field = pick(field$gp_field_si),
                 bathy_effect = pick(field$bathy_effect_si),
                 spline_setup = field$bathy_spline, gp_params = gp_params[lt_species],
                 beta_bathy_true = field$beta_bathy_true[match(lt_species, names(gp_params)), , drop = FALSE],
                 sp_params = sp_params[lt_species]),
    observed = observed_lt)

  list(edna = edna, visual = visual,
       field = c(field, list(locations = locs)),
       meta = list(seed = seed, domain = dom, lt_species = lt_species))
}

# Sanity check that the two sources share one field: for station/segment
# pairs closer than `max_km`, the true log density (per species) should agree
# up to the GP's short-range variation and the independent bottom-depth draws.
joint_field_agreement <- function(sim, max_km = 10) {
  L <- sim$field$locations
  st <- which(L$source == "station"); sg <- which(L$source == "segment")
  d <- sqrt(outer(L$X[st], L$X[sg], "-")^2 + outer(L$Y[st], L$Y[sg], "-")^2)
  near <- which(d < max_km, arr.ind = TRUE)
  if (!nrow(near)) return(NULL)
  i <- st[near[, 1]]; j <- sg[near[, 2]]
  f <- sim$field$gp_field_si
  data.frame(species = sim$field$species, n_pairs = nrow(near),
             cor_gp_field = vapply(seq_len(ncol(f)), function(s) cor(f[i, s], f[j, s]), numeric(1)))
}

# -----------------------------------------------------------------------------
# format_stan_data_joint(): Stan data for stan/hsgp_joint.stan from a
# simulate_joint() object (or any object with the same $edna / $visual shape).
#
#   use_edna, use_visual  which likelihoods enter (1/1 joint; 1/0 eDNA-only;
#                         0/1 visual-only). Everything else - field, priors,
#                         basis - is identical across the three, so they differ
#                         only in the data used.
#   conv_free             species names whose eDNA conversion factor is
#                         estimated (log_conv ~ normal(prior_log_conv)); the
#                         rest stay fixed at sim$edna$meta$conv_factor.
#   prior_mu_sp           list(mu, sig) for the log-density intercepts, by
#                         species name (+ "junk"). Default: the visual model's
#                         normal(-5, 2) for the line-transect species and
#                         HSGP4eDNA's normal(2, 1.5) for the others and junk.
#
# Coordinates of both sources are normalised by the DOMAIN extents, and the
# design-based GP prior is built from the union of station and segment
# locations (sample locations only - no detections).
# The eDNA observation data come from HSGP4eDNA's format_stan_data_gp2d_bathysp()
# and the visual data from format_stan_data_visual() (single ESW: no sea state;
# its detection key, hazard-rate prior and per-species truncation carry over).
# -----------------------------------------------------------------------------
format_stan_data_joint <- function(sim, use_edna = 1L, use_visual = 1L,
                                   conv_free = character(0),
                                   prior_log_conv = c(mu = log(100), sig = 2),
                                   prior_mu_sp = NULL,
                                   HSGP_M = NULL, HSGP_C = c(1.5, 1.5),
                                   M_max = GP_BASIS_BUDGET) {
  dom    <- sim$meta$domain
  lt     <- sim$meta$lt_species
  sp_all <- names(sim$edna$meta$sp_names)            # hake, humpback, pwsd
  S      <- length(sp_all)
  coord_centre <- c(dom$X_km_max, dom$Y_km_max) / 2
  coord_scale  <- coord_centre
  norm <- function(xy) sweep(sweep(as.matrix(xy), 2, coord_centre, "-"), 2, coord_scale, "/")

  # eDNA observation data (its own coords / GP prior / basis are replaced below)
  ed  <- format_stan_data_gp2d_bathysp(sim$edna, prediction_n = 2L)$stan_data
  smp <- sim$edna$design$samples
  seg <- sim$visual$design$segments

  # Design-based GP prior from the union of station and segment locations
  design_xy <- unique(rbind(as.matrix(smp[, c("X", "Y")]), as.matrix(seg[, c("X", "Y")])))
  gp_prior <- gp_prior_from_design(design_xy, domain_range = c(dom$X_km_max, dom$Y_km_max),
                                   coord_scale = coord_scale, M_max = M_max, c = HSGP_C[1])
  if (is.null(HSGP_M)) HSGP_M <- gp_prior$HSGP_M
  INDICES <- as.matrix(do.call(tidyr::expand_grid, lapply(HSGP_M, seq_len)))

  # Visual data per line-transect species (single ESW). The sea-state
  # covariate is not in the joint model yet (roadmap Phase 6).
  wants_bft <- vapply(lt, function(sp) isTRUE(sim$visual$truth$sp_params[[sp]]$use_bft), logical(1))
  if (any(wants_bft))
    warning("the joint model has no sea-state covariate yet; ignoring use_bft for: ",
            paste(lt[wants_bft], collapse = ", "), call. = FALSE)
  vd <- lapply(setNames(lt, lt), function(sp)
    format_stan_data_visual(sim$visual, sp, HSGP_M = HSGP_M, use_bft = FALSE))
  stopifnot(isTRUE(all.equal(vd[[1]]$B_bathy, vd[[length(vd)]]$B_bathy)),
            ncol(vd[[1]]$B_bathy) == ed$K_bathy)
  pick <- function(field) vapply(vd, function(v) as.numeric(v[[field]]), numeric(1))
  picki <- function(field) vapply(vd, function(v) as.integer(v[[field]]), integer(1))

  # Intercept priors (species + junk)
  if (is.null(prior_mu_sp)) {
    mu  <- setNames(rep(2.0, S + 1), c(sp_all, "junk"))
    sig <- setNames(rep(1.5, S + 1), c(sp_all, "junk"))
    mu[lt] <- -5; sig[lt] <- 2
    prior_mu_sp <- list(mu = mu, sig = sig)
  }

  conv_idx <- match(conv_free, sp_all)
  stopifnot(!anyNA(conv_idx))
  noise_on <- as.integer(use_edna == 1L)

  sd <- list(
    use_edna = as.integer(use_edna), use_visual = as.integer(use_visual),
    # field
    S = S, D1 = 2L, M = as.integer(prod(HSGP_M)), INDICES = INDICES, L_hsgp = HSGP_C,
    coord_scale = coord_scale, K_bathy = ed$K_bathy,
    prior_gp_sigma_shape = gp_prior$prior_gp_sigma_shape,
    prior_gp_sigma_rate  = gp_prior$prior_gp_sigma_rate,
    prior_gp_l_shape = gp_prior$prior_gp_l_shape, prior_gp_l_scale = gp_prior$prior_gp_l_scale,
    prior_mu_sp_mu = as.numeric(prior_mu_sp$mu[c(sp_all, "junk")]),
    prior_mu_sp_sig = as.numeric(prior_mu_sp$sig[c(sp_all, "junk")]),
    prior_beta_bathy_sig = ed$prior_beta_bathy_sig,
    # eDNA
    N = ed$N, coords = norm(smp[, c("X", "Y")]), B_bathy = ed$B_bathy,
    log_zsample_effect = ed$log_zsample_effect, log_conv_factor = ed$log_conv_factor,
    log_vol_filtered = ed$log_vol_filtered, vol_aliquot = ed$vol_aliquot,
    N_qpcr_long = ed$N_qpcr_long, qpcr_sample_idx = ed$qpcr_sample_idx,
    qpcr_detect = ed$qpcr_detect, qpcr_ct = ed$qpcr_ct,
    N_mb_long = ed$N_mb_long, mb_sample_idx = ed$mb_sample_idx,
    mb_reads = ed$mb_reads, mb_total = ed$mb_total,
    alpha_ct = ed$alpha_ct, beta_ct = ed$beta_ct, kappa = ed$kappa,
    gamma0_ct = ed$gamma0_ct, gamma1_ct = ed$gamma1_ct, sigma0_ct = ed$sigma0_ct,
    prior_beta0_phi_mu = ed$prior_beta0_phi_mu, prior_beta0_phi_sig = ed$prior_beta0_phi_sig,
    prior_gamma0_phi_mu = ed$prior_gamma0_phi_mu, prior_gamma0_phi_sig = ed$prior_gamma0_phi_sig,
    prior_gamma1_phi_mu = ed$prior_gamma1_phi_mu, prior_gamma1_phi_sig = ed$prior_gamma1_phi_sig,
    # noise terms only when the eDNA data are used (else they would be
    # prior-only latent variables)
    use_bottle_noise = ed$use_bottle_noise * noise_on,
    bottle_noise_dist = ed$bottle_noise_dist, bottle_theta_shared = ed$bottle_theta_shared,
    prior_cv_bottle_sig = ed$prior_cv_bottle_sig,
    use_junk_noise = ed$use_junk_noise * noise_on, prior_cv_junk_sig = ed$prior_cv_junk_sig,
    n_conv_free = length(conv_idx), conv_free_idx = as.array(as.integer(conv_idx)),
    prior_log_conv_mu = as.array(rep(prior_log_conv[["mu"]], length(conv_idx))),
    prior_log_conv_sig = as.array(rep(prior_log_conv[["sig"]], length(conv_idx))),
    # visual
    L = length(lt), lt_species_idx = as.array(match(lt, sp_all)),
    n_seg = nrow(seg), seg_coords = norm(seg[, c("X", "Y")]), B_bathy_seg = vd[[1]]$B_bathy,
    seg_l = as.numeric(seg$seg_l),
    seg_count = do.call(rbind, lapply(vd, function(v) as.integer(v$seg_count))),
    w = as.array(pick("w")),
    n_det = as.array(picki("n")), n_det_total = sum(picki("n")),
    det_x = as.array(unlist(lapply(vd, `[[`, "x"), use.names = FALSE)),
    det_s = as.array(as.integer(unlist(lapply(vd, `[[`, "s"), use.names = FALSE))),
    log_sigma_prior_mean = as.array(pick("log_sigma_prior_mean")),
    log_sigma_prior_sd = as.array(pick("log_sigma_prior_sd")),
    use_size_covar = as.array(picki("use_size_covar")), s_centre = as.array(pick("s_centre")),
    beta_size_prior_mean = as.array(pick("beta_size_prior_mean")),
    beta_size_prior_sd = as.array(pick("beta_size_prior_sd")),
    det_key = as.array(picki("det_key")),
    log_b_prior_mean = as.array(pick("log_b_prior_mean")),
    log_b_prior_sd = as.array(pick("log_b_prior_sd")),
    S_max = as.array(picki("S_max")), group_size_dist = as.array(picki("group_size_dist")),
    mu_s_prior_shape = as.array(pick("mu_s_prior_shape")),
    mu_s_prior_rate = as.array(pick("mu_s_prior_rate")),
    phi_s_prior_shape = as.array(pick("phi_s_prior_shape")),
    phi_s_prior_rate = as.array(pick("phi_s_prior_rate")),
    mu_log_prior_mean = as.array(pick("mu_log_prior_mean")),
    mu_log_prior_sd = as.array(pick("mu_log_prior_sd")),
    sigma_log_prior_shape = as.array(pick("sigma_log_prior_shape")),
    sigma_log_prior_rate = as.array(pick("sigma_log_prior_rate"))
  )
  attr(sd, "gp_prior") <- gp_prior
  attr(sd, "species") <- sp_all
  sd
}
