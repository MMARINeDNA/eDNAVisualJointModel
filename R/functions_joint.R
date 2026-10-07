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
                           n_transects = 25L, seg_length = 10, w = 5.0, df = 4L,
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
