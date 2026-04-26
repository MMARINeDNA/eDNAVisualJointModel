# =============================================================================
# simulate_whale_edna_v3.2.R
#
# Single-species (Pacific hake), surface-only (Z_sample = 0), qPCR-only
# debug case for the v3 HSGP pipeline.
#
# Domain, bathymetry, station design, and hake habitat preference are
# identical to v3. This script strips the model down to its simplest
# form so we can diagnose whether the HSGP latent field is recoverable
# at all from a single observation type.
#
# Changes from v3:
#   * n_species = 1 (hake only). Humpback / PWSD removed.
#   * sample_depths = 0 only. Single sample per station at the surface;
#     Z_sample is always 0 so the water-column effect is the constant
#     exp(0) = 1 (no contribution).
#   * Metabarcoding generation removed. qPCR observations only.
#   * Output is otherwise structurally identical to v3 (same `sim`
#     list shape, just with S=1 slices), so downstream scripts can
#     treat it as a degenerate v3 sim with the MB block ignored.
# =============================================================================

library(tidyverse)
library(MASS)
set.seed(42)

# ---------------------------------------------------------------------------
# 0. Study design - hake only, surface only
# ---------------------------------------------------------------------------
n_species  <- 1
sp_names   <- c("Merluccius_productus")
sp_common  <- c("Pacific hake")

has_qpcr   <- c(TRUE)
n_qpcr_rep <- 3
R <- n_qpcr_rep
marker     <- "MARVER1"

# UTM Zone 10N bounding box (metres): San Francisco (~37.77 deg N) to 49 deg N
X_min <- 100000;  X_max <- 600000   # 500 km E-W
Y_min <- 4180000; Y_max <- 5450000  # ~1270 km N-S

X_km_max <- (X_max - X_min) / 1000   # 500 km
Y_km_max <- (Y_max - Y_min) / 1000   # 1270 km

n_stations    <- 300
conv_factor   <- 10
sample_depths <- c(0)               # surface only
n_sample_depth <- length(sample_depths)

vol_filtered  <- 2.5
vol_aliquot   <- 2

# ---------------------------------------------------------------------------
# 1. Rotated bathymetry (unchanged from v3)
# ---------------------------------------------------------------------------
rotation_deg <- 25
rot <- rotation_deg * pi / 180

X_prime_fn <- function(X_km, Y_km) {
  raw <- cos(rot) * (X_km / X_km_max) + sin(rot) * (Y_km / Y_km_max)
  300 * raw / (cos(rot) + sin(rot))
}

bathy_mean_fn <- function(X_km, Y_km) {
  xp       <- X_prime_fn(X_km, Y_km)
  abyssal  <- 2500
  shelf    <- 80
  slope    <- abyssal + (shelf - abyssal) / (1 + exp(-0.06 * (xp - 180)))
  pmax(slope, 10)
}

# ---------------------------------------------------------------------------
# 2. Station locations (unchanged from v3)
# ---------------------------------------------------------------------------
set.seed(1)

gen_stratum <- function(n, xp_lo, xp_hi) {
  out <- data.frame(X_km = numeric(0), Y_km = numeric(0))
  while (nrow(out) < n) {
    cand <- data.frame(
      X_km = runif(n * 4, 0, X_km_max),
      Y_km = runif(n * 4, 0, Y_km_max)
    )
    cand$xp <- X_prime_fn(cand$X_km, cand$Y_km)
    keep <- cand[cand$xp >= xp_lo & cand$xp <= xp_hi, c("X_km", "Y_km")]
    out <- rbind(out, keep)
  }
  out[seq_len(n), ]
}

n_shelf    <- round(n_stations * 0.50)
n_slope    <- round(n_stations * 0.30)
n_offshore <- n_stations - n_shelf - n_slope

stations <- bind_rows(
  gen_stratum(n_shelf,    xp_lo = 220, xp_hi = 300),
  gen_stratum(n_slope,    xp_lo = 120, xp_hi = 220),
  gen_stratum(n_offshore, xp_lo =   0, xp_hi = 120)
) %>%
  slice_sample(n = n_stations) %>%
  mutate(
    station = row_number(),
    X       = X_km,
    Y       = Y_km,
    X_utm   = X_min + X_km * 1000,
    Y_utm   = Y_min + Y_km * 1000,
    X_prime = X_prime_fn(X, Y)
  ) %>%
  dplyr::select(-X_km, -Y_km)

set.seed(2)
stations <- stations %>%
  mutate(
    Z_bathy_mean  = bathy_mean_fn(X, Y),
    Z_bathy_noise = rnorm(n_stations, 0, Z_bathy_mean * 0.12),
    Z_bathy       = pmax(Z_bathy_mean + Z_bathy_noise, 10)
  )

cat(sprintf("Stations: %d\n", n_stations))
cat(sprintf("  X_prime range: %.0f - %.0f km\n",
            min(stations$X_prime), max(stations$X_prime)))
cat(sprintf("  Z_bathy range: %.0f - %.0f m\n",
            min(stations$Z_bathy), max(stations$Z_bathy)))

# Surface-only sample table (one row per station)
samples <- expand.grid(
  station   = 1:n_stations,
  depth_idx = 1:n_sample_depth
) %>%
  left_join(stations, by = "station") %>%
  mutate(
    Z_sample  = sample_depths[depth_idx],
    sample_id = row_number()
  ) %>%
  filter(Z_sample <= Z_bathy) %>%
  mutate(sample_id = row_number())

N <- nrow(samples)
cat(sprintf("Samples (surface only): %d\n", N))

coords_gp <- as.matrix(samples[, c("X", "Y", "Z_bathy")])

# ---------------------------------------------------------------------------
# 3. GP hyperparameters - hake only
# ---------------------------------------------------------------------------

gp_params <- list(

  hake = list(
    name  = sp_names[1],
    sigma = 1.2,
    lx    =  50,
    ly    = 300,
    lz    = 150,
    mu    = log(20),
    # Surface only: pref at z=0 is 0.0, so exp(0) = 1 (multiplier = 1)
    zsample_pref = c(0.0)
    # Habitat-preference fields (zbathy_pref_*, y_pref_*) removed in v3.2
    # Option A: the simulation now matches the model's structural
    # assumption that f_s is a zero-mean GP. The deterministic preference
    # mean function present in v3 was being absorbed into the model's GP
    # field, making length-scale recovery against the named (lx, ly, lz)
    # ambiguous because the model was estimating the length-scale of the
    # combined preference + residual structure.
  )

)

# ---------------------------------------------------------------------------
# 4. Preference + GP draws
# ---------------------------------------------------------------------------

aniso_cov <- function(coords, sigma, lx, ly, lz) {
  n <- nrow(coords)
  K <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in i:n) {
      d2 <- ((coords[i,1]-coords[j,1])/lx)^2 +
        ((coords[i,2]-coords[j,2])/ly)^2 +
        ((coords[i,3]-coords[j,3])/lz)^2
      K[i,j] <- K[j,i] <- sigma^2 * exp(-0.5 * d2)
    }
  }
  K + diag(1e-6, n)
}

cat("Drawing zero-mean GP realisation (N =", N, ")...\n")
lambda_true_si <- matrix(NA, N, n_species)
gp_field_si    <- matrix(NA, N, n_species)
gp_mean_si     <- matrix(0,  N, n_species)   # zero-mean GP - kept for output shape

for (s in seq_len(n_species)) {
  p  <- gp_params[[s]]

  K  <- aniso_cov(coords_gp, p$sigma, p$lx, p$ly, p$lz)

  fs <- as.vector(mvrnorm(1, mu = rep(0.0, N), Sigma = K))
  gp_field_si[, s] <- fs

  lambda_true_si[, s] <- exp(p$mu + fs)
}

# ---------------------------------------------------------------------------
# 5. eDNA water column process (trivial at surface)
# ---------------------------------------------------------------------------

zsample_effect <- matrix(1.0, N, n_species)   # exp(0) at surface

C_obs_si        <- matrix(NA, N, n_species)
C_obs_si_a_qpcr <- matrix(NA, N*R, n_species + 1)
C_obs_si_a_qpcr[, n_species + 1] <- rep(1:N, each = R)

for (s in seq_len(n_species)){
  C_obs_si[, s] <- MASS::rnegbin(
    N,
    mu    = conv_factor * lambda_true_si[, s] * zsample_effect[, s] * vol_filtered,
    theta = 10
  )
  C_obs_si_a_qpcr[, s] <- rbinom(N*R, rep(C_obs_si[, s], each = R), vol_aliquot / 100)
}

# ---------------------------------------------------------------------------
# 6. qPCR hurdle - hake (only species)
# ---------------------------------------------------------------------------
qpcr_p <- list(
  kappa    = 0.85,
  alpha_ct = 38.0,
  beta_ct  = 1.44,
  sigma_ct = 0.50
)

qpcr_detect_raw <- matrix(NA_integer_, N * R, n_species)
qpcr_ct_raw     <- matrix(NA_real_,    N * R, n_species)

p_det_hake <- 1 - exp(-qpcr_p$kappa * C_obs_si_a_qpcr[, 1])
det <- rbinom(N * R, 1, p_det_hake)
Ct  <- ifelse(
  det == 1L,
  qpcr_p$alpha_ct -
    qpcr_p$beta_ct * log(pmax(C_obs_si_a_qpcr[, 1], 1)) +
    rnorm(N * R, 0, qpcr_p$sigma_ct),
  NA_real_
)
qpcr_detect_raw[, 1] <- det
qpcr_ct_raw[, 1]     <- Ct

cat(sprintf("qPCR detection rate (hake): %.1f%%\n", 100 * mean(det)))

# ---------------------------------------------------------------------------
# 7. Bundle and save
# ---------------------------------------------------------------------------
sim <- list(
  meta = list(
    n_species      = n_species,
    sp_names       = sp_names,
    sp_common      = sp_common,
    has_qpcr       = has_qpcr,
    n_qpcr_rep     = n_qpcr_rep,
    marker         = marker,
    n_stations     = n_stations,
    sample_depths  = sample_depths,
    n_sample_depth = n_sample_depth,
    vol_filtered   = vol_filtered,
    N              = N,
    utm_zone       = "10N",
    X_min_utm      = X_min,
    X_max_utm      = X_max,
    Y_min_utm      = Y_min,
    Y_max_utm      = Y_max,
    X_km_max       = X_km_max,
    Y_km_max       = Y_km_max,
    rotation_deg   = rotation_deg,
    bathy_profile  = list(
      rotation_deg       = rotation_deg,
      shelf_break_xprime = 180,
      abyssal_depth      = 2500,
      shelf_depth        =   80,
      slope_k            = 0.06
    ),
    coord_note = paste(
      "X (km): 0 = western edge of domain (100000 UTM E),",
      "500 = eastern edge (600000 UTM E); ",
      "Y (km): 0 = SF (~4180000 UTM N, 37.77 deg N),",
      "1270 = US/Canada border (~5450000 UTM N, ~49 deg N); ",
      "Z_bathy: bottom depth (m), derived from rotated X_prime; ",
      "Z_sample: surface only (0 m)"
    )
  ),
  design = list(
    stations = stations,
    samples  = samples
  ),
  truth = list(
    gp_field_si      = gp_field_si,
    gp_mean_si       = gp_mean_si,
    lambda_true_si   = lambda_true_si,
    C_obs_si         = C_obs_si,
    C_obs_si_a_qpcr  = C_obs_si_a_qpcr,
    zsample_effect   = zsample_effect,
    gp_params        = gp_params,
    qpcr_params      = qpcr_p
  ),
  observed = list(
    qpcr_detect_raw = qpcr_detect_raw,
    qpcr_ct_raw     = qpcr_ct_raw
  )
)

dir.create("outputs/whale_edna_output_v3.2", showWarnings = FALSE, recursive = TRUE)
saveRDS(sim, "outputs/whale_edna_output_v3.2/whale_edna_sim_v3.2.rds")
cat("Saved outputs/whale_edna_output_v3.2/whale_edna_sim_v3.2.rds\n")
