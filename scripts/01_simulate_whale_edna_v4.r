# =============================================================================
# simulate_whale_edna_v4.R
#
# Anisotropic 3-species eDNA simulation
# Domain: US West Coast, San Francisco to US/Canada border, UTM Zone 10N
#
#   X (Easting)  : 100,000 – 600,000 m   (500 km)
#   Y (Northing) : 4,180,000 – 5,450,000 m (~1270 km, SF to 49°N)
#   Z_bathy      : bottom depth at each station (m), derived from a ROTATED
#                  cross-shore coordinate.
#
# Changes from v3:
#   * **Variable per-sample replication**. Most samples have the full 3
#     qPCR + 3 metabarcoding replicates, but ~20% of samples (chosen
#     independently for each process) are reduced to 1 or 2 reps,
#     matching the loss / QC failure rate seen in real surveys.
#   * qPCR and MB observed data are now stored in **long form** directly
#     (one row per sample-aliquot), with explicit `qpcr_sample_idx` /
#     `mb_sample_idx` vectors and per-sample rep counts. The format
#     script no longer has to re-derive the long-form structure.
#   * **Junk read category** on the MB side — a lumped background of
#     non-target species / primer artifacts that claims a random 5–95%
#     share of each MB aliquot's read depth. Target-species detection
#     now reflects this competition realistically. Junk does not appear
#     on the qPCR side (qPCR primers are species-specific).
#   * **Species densities retuned** so detection rates match real-data
#     reference values with MARVER1: hake ≈ 95% in both qPCR and MB,
#     humpback ≈ 20% in MB, PWS dolphin ≈ 20% in MB (with independent
#     spatial distributions — they aren't necessarily the same 20%).
#
# Two distinct depth processes (unchanged from v3):
#   1. Z_bathy  — bathymetric depth at (X, Y); drives TRUE animal density
#                  via species habitat preferences.
#   2. Z_sample — depth of the water sample; drives eDNA vertical
#                  distribution within the water column.
# =============================================================================

library(tidyverse)
library(MASS)
set.seed(42)

# ---------------------------------------------------------------------------
# 0. Study design
# ---------------------------------------------------------------------------
n_species  <- 3
sp_names   <- c("Merluccius_productus",
                "Megaptera_novaeangliae",
                "Lagenorhynchus_obliquidens")
sp_common  <- c("Pacific hake",
                "Humpback whale",
                "Pacific white-sided dolphin")

has_qpcr   <- c(TRUE,  FALSE, FALSE)
n_qpcr_rep <- 3            # full-design qPCR replicates per sample
n_mb_rep   <- 3            # full-design metabarcoding replicates per sample
marker     <- "MARVER1"

# Fraction of samples where a process is reduced below its full 3 reps.
# A "reduced" sample keeps 1 or 2 reps (chosen uniformly). Applied
# independently to qPCR and MB, so a sample can be reduced in one process
# and not the other.
prop_reduced_qpcr_rep <- 0.20
prop_reduced_mb_rep   <- 0.20

# UTM Zone 10N bounding box (metres): San Francisco (~37.77°N) to 49°N
X_min <- 100000;  X_max <- 600000   # 500 km E-W
Y_min <- 4180000; Y_max <- 5450000  # ~1270 km N-S

# Domain extent in km (for all GP / preference calculations)
X_km_max <- (X_max - X_min) / 1000   # 500 km
Y_km_max <- (Y_max - Y_min) / 1000   # 1270 km

n_stations    <- 300
conv_factor <- 10
# Three sample depths per station: surface, mid-column, upper slope
sample_depths <- c(0, 150, 500)
n_sample_depth <- length(sample_depths)

vol_filtered  <- 2.5
vol_aliquot <- 2
n_mb_reads    <- 50000

# ---------------------------------------------------------------------------
# 1. Rotated bathymetry
#
# Rather than depth = f(X), use a rotated cross-shore axis X_prime so the
# shelf-slope isobaths run NW-SE (roughly parallel to the US West Coast).
#
# The rotation is applied in NORMALIZED space (both axes rescaled to [0, 1])
# so the angle has a consistent geometric meaning regardless of the domain
# aspect ratio. At 0° the bathymetry depends only on X (V2 behaviour); at
# 90° only on Y. We use 25° — large enough for clearly tilted isobaths, but
# small enough that both shelf and offshore are present at every latitude.
# (At 45° in this elongated domain, latitude and bathymetry become almost
# perfectly confounded, which breaks realistic species preferences.)
# ---------------------------------------------------------------------------
rotation_deg <- 25
rot <- rotation_deg * pi / 180

# X_prime spans [0, 300] km-equivalent over the domain, matching the V2
# shelf profile. 0 = most offshore (SW corner), 300 = most nearshore (NE corner).
X_prime_fn <- function(X_km, Y_km) {
  raw <- cos(rot) * (X_km / X_km_max) + sin(rot) * (Y_km / Y_km_max)
  300 * raw / (cos(rot) + sin(rot))
}

# Shelf-slope profile on the rotated axis — same logistic shape as V2,
# now driven by X_prime instead of raw X
bathy_mean_fn <- function(X_km, Y_km) {
  xp       <- X_prime_fn(X_km, Y_km)
  abyssal  <- 2500
  shelf    <- 80
  slope    <- abyssal + (shelf - abyssal) / (1 + exp(-0.06 * (xp - 180)))
  pmax(slope, 10)
}

# ---------------------------------------------------------------------------
# 2. Station locations
#
# Stratified by the rotated cross-shore coordinate so we have dense
# coverage of the shelf/slope and lighter coverage offshore.
# ---------------------------------------------------------------------------
set.seed(1)

gen_stratum <- function(n, xp_lo, xp_hi) {
  # Rejection sampling: draw uniformly in (X_km, Y_km) and keep points
  # whose rotated X_prime falls in [xp_lo, xp_hi]
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

# Three strata matched to the shelf profile:
#   shelf    : X_prime > 220 km (Z_bathy ~  80 – 300 m)
#   slope    : X_prime 120–220 (Z_bathy ~ 300 – 2000 m)
#   offshore : X_prime <  120  (Z_bathy > 2000 m)
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
cat(sprintf("  X_prime range: %.0f – %.0f km\n",
            min(stations$X_prime), max(stations$X_prime)))
cat(sprintf("  Z_bathy range: %.0f – %.0f m\n",
            min(stations$Z_bathy), max(stations$Z_bathy)))

# Full sample table: station × sample depth, drop samples below seafloor
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
cat(sprintf("Samples after filtering Z_sample > Z_bathy: %d\n", N))

# GP coordinates: (X, Y, Z_bathy). Bathymetry drives density structure.
coords_gp <- as.matrix(samples[, c("X", "Y", "Z_bathy")])

# ---------------------------------------------------------------------------
# 3. GP hyperparameters + realistic spatial structure
#
# Model form exposed to downstream inference:
#   log(lambda_s) = mu_s + f_s(X, Y, Z_bathy)
#
# All spatial structure (bathymetric habitat, latitude band) lives inside
# f_s — the GP. In the simulation, f_s is drawn from an anisotropic
# Gaussian with a non-zero mean that encodes the species' habitat
# preferences, and an anisotropic squared-exponential covariance. From the
# model's perspective this is one realisation of a zero-mean GP; the
# habitat structure is simply the feature we want the data to reveal.
#
# Realistic preferences for the US West Coast:
#   Hake         — shelf-slope break (Z_bathy 150–400 m), broad latitude
#   Humpback     — shallow shelf, concentrated south (central CA to S. OR)
#   PWS dolphin  — deep slope / oceanic, concentrated north of Cape Mendocino
# ---------------------------------------------------------------------------

gp_params <- list(

  hake = list(
    name  = sp_names[1],
    sigma = 1.2,
    lx    =  50,
    ly    = 300,
    lz    = 150,
    # mu chosen so hake appears in ~95% of qPCR and MB samples (real-data
    # target). Hake density is much higher than the whales below.
    mu    = log(6),
    # Bathymetric preference: shelf-slope break
    zbathy_pref_mu  = 200,
    zbathy_pref_sd  = 100,
    zbathy_pref_amp = 2.0,
    # Latitude preference: broad, slightly centred on OR-WA shelf
    y_pref_mu  = 800,      # km north of SF (~central OR)
    y_pref_sd  = 500,
    y_pref_amp = 1.0,
    # eDNA water column at (0, 150, 500 m): peaks at 150 m (mid-column)
    zsample_pref = c( 0.0,  3.0,  1.5)
  ),

  humpback = list(
    name  = sp_names[2],
    sigma = 1.0,
    lx    =  50,
    ly    = 300,
    lz    = 100,
    # mu chosen so humpback appears in ~20% of MB samples (real-data
    # target with MARVER1). Independent from PWSD. Two-plus orders of
    # magnitude rarer than hake, reflecting real animal density as well
    # as lower eDNA share when hake dominates the multinomial.
    mu    = log(0.05),
    # Bathymetric preference: nearshore shelf (50–150 m)
    zbathy_pref_mu  = 100,
    zbathy_pref_sd  =  80,
    zbathy_pref_amp = 2.0,
    # Latitude preference: SOUTHERN (SF / central CA / S. OR)
    y_pref_mu  = 150,      # km north of SF
    y_pref_sd  = 350,
    y_pref_amp = 2.0,
    # eDNA water column at (0, 150, 500 m): peaks at feeding-dive depths
    zsample_pref = c( 0.0,  1.1, -1.3)
  ),

  pwsd = list(
    name  = sp_names[3],
    sigma = 1.3,
    lx    =  40,
    ly    = 300,
    lz    = 300,
    # mu chosen so PWSD appears in ~20% of MB samples (real-data target
    # with MARVER1). Independent spatial distribution from humpback.
    mu    = log(0.12),
    # Bathymetric preference: offshore slope / oceanic (>500 m)
    zbathy_pref_mu  = 800,
    zbathy_pref_sd  = 400,
    zbathy_pref_amp = 2.0,
    # Latitude preference: NORTHERN (N. OR / WA / Vancouver I. shelf edge)
    y_pref_mu  = 1050,     # km north of SF (~N. WA)
    y_pref_sd  = 350,
    y_pref_amp = 2.0,
    # eDNA water column at (0, 150, 500 m): surface-active, drops off with depth
    zsample_pref = c( 0.0, -2.0, -3.5)
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

# Gaussian-shaped preference over a covariate, scaled to [-amp/2, +amp/2]
gauss_pref <- function(x, mu, sd, amp) {
  raw <- dnorm(x, mu, sd)
  raw <- raw / dnorm(mu, mu, sd)   # peak = 1
  amp * (raw - 0.5)
}

zbathy_pref_fn <- function(Z_bathy, mu_z, sd_z, amp)
  gauss_pref(Z_bathy, mu_z, sd_z, amp)

y_pref_fn <- function(Y, mu_y, sd_y, amp)
  gauss_pref(Y, mu_y, sd_y, amp)

cat("Drawing GP realisations (N =", N, ")...\n")
lambda_true_si       <- matrix(NA, N, n_species)
gp_field_si          <- matrix(NA, N, n_species)
gp_mean_si           <- matrix(NA, N, n_species)   # diagnostic: GP mean only

for (s in seq_len(n_species)) {
  p  <- gp_params[[s]]

  K  <- aniso_cov(coords_gp, p$sigma, p$lx, p$ly, p$lz)

  # Non-zero GP mean encoding species habitat preferences. The GP draw
  # below folds this into f_s, so the exposed model is log(lambda) = mu + f.
  gp_mean <- zbathy_pref_fn(samples$Z_bathy,
                            p$zbathy_pref_mu,
                            p$zbathy_pref_sd,
                            p$zbathy_pref_amp) +
             y_pref_fn(samples$Y,
                       p$y_pref_mu,
                       p$y_pref_sd,
                       p$y_pref_amp)
  gp_mean_si[, s] <- gp_mean

  fs <- as.vector(mvrnorm(1, mu = gp_mean, Sigma = K))
  gp_field_si[, s] <- fs

  lambda_true_si[, s] <- exp(p$mu + fs)
}

# ---------------------------------------------------------------------------
# 5. eDNA water column process
# ---------------------------------------------------------------------------

zsample_effect <- matrix(NA, N, n_species)
for (s in seq_len(n_species)) {
  p <- gp_params[[s]]
  pref_fn <- approxfun(
    x      = sample_depths,
    y      = p$zsample_pref,
    method = "linear",
    rule   = 2
  )
  zsample_effect[, s] <- exp(pref_fn(samples$Z_sample))
}

C_obs_si <- matrix(NA_integer_, N, n_species)
for (s in 1:n_species) {
  C_obs_si[, s] <- MASS::rnegbin(
    N,
    mu    = conv_factor * lambda_true_si[, s] * zsample_effect[, s] * vol_filtered,
    theta = 10
  )
}

# ---------------------------------------------------------------------------
# 6. Per-sample replicate counts (variable replication)
# ---------------------------------------------------------------------------
# For each sample, each process is independently either at full design
# (n_{qpcr,mb}_rep) or reduced to {1, 2} reps with probability
# prop_reduced_{qpcr,mb}_rep.
set.seed(99)
reduced_qpcr  <- runif(N) < prop_reduced_qpcr_rep
n_qpcr_rep_i  <- ifelse(reduced_qpcr,
                        sample.int(2L, N, replace = TRUE),   # {1, 2}
                        as.integer(n_qpcr_rep))
reduced_mb    <- runif(N) < prop_reduced_mb_rep
n_mb_rep_i    <- ifelse(reduced_mb,
                        sample.int(2L, N, replace = TRUE),
                        as.integer(n_mb_rep))

N_qpcr_long <- sum(n_qpcr_rep_i)
N_mb_long   <- sum(n_mb_rep_i)

qpcr_sample_idx <- rep(seq_len(N), times = n_qpcr_rep_i)
mb_sample_idx   <- rep(seq_len(N), times = n_mb_rep_i)

cat(sprintf("qPCR reps: full=%d reduced=%d (%.0f%%)  long rows=%d\n",
            sum(!reduced_qpcr), sum(reduced_qpcr),
            100 * mean(reduced_qpcr), N_qpcr_long))
cat(sprintf("MB   reps: full=%d reduced=%d (%.0f%%)  long rows=%d\n",
            sum(!reduced_mb), sum(reduced_mb),
            100 * mean(reduced_mb),   N_mb_long))

# ---------------------------------------------------------------------------
# 7. Aliquot copies (long form)
# ---------------------------------------------------------------------------
# Each qPCR / MB aliquot is a Binomial draw from the bottle copy count
# C_obs_si, with probability vol_aliquot / 100 (2 µL aliquot from a 100 µL
# elution).
qpcr_copies <- rbinom(N_qpcr_long,
                      size = C_obs_si[qpcr_sample_idx, 1],   # hake only
                      prob = vol_aliquot / 100)

mb_copies <- matrix(NA_integer_, N_mb_long, n_species)
for (s in seq_len(n_species)) {
  mb_copies[, s] <- rbinom(N_mb_long,
                           size = C_obs_si[mb_sample_idx, s],
                           prob = vol_aliquot / 100)
}

# ---------------------------------------------------------------------------
# 8. qPCR hurdle — hake only
# ---------------------------------------------------------------------------
qpcr_p <- list(
  kappa    = 0.85,
  alpha_ct = 38.0,
  beta_ct  = 1.44,
  sigma_ct = 0.50
)

p_det_hake   <- 1 - exp(-qpcr_p$kappa * qpcr_copies)
qpcr_detect  <- rbinom(N_qpcr_long, 1, p_det_hake)
qpcr_ct      <- ifelse(
  qpcr_detect == 1L,
  qpcr_p$alpha_ct -
    qpcr_p$beta_ct * log(pmax(qpcr_copies, 1)) +
    rnorm(N_qpcr_long, 0, qpcr_p$sigma_ct),
  NA_real_
)

# ---------------------------------------------------------------------------
# 9. Metabarcoding — "junk" background + target-species multinomial
# ---------------------------------------------------------------------------
# Real MB data is dominated by non-target species (plankton, bacteria,
# other vertebrates that MARVER1 picks up, primer artifacts, etc.). We
# model this with a single lumped "junk" category that:
#   * is present in every MB aliquot (but not in qPCR)
#   * claims a random 5–95% share of the reads in each aliquot
#   * competes with the target species for read budget, so the three
#     target species' reads are drawn from the remaining (100% – junk%)
#     of the depth.
# The Stan model still sees only the three target species: mb_total is
# the sum of target reads per aliquot, and pi_edna is the target-species
# proportion.
set.seed(200)
read_depth <- round(runif(N_mb_long, n_mb_reads/2, n_mb_reads))
pi_junk    <- runif(N_mb_long, 0.05, 0.95)
junk_reads <- rbinom(N_mb_long, size = read_depth, prob = pi_junk)
target_read_depth <- read_depth - junk_reads

# Target-species proportions (within the non-junk pool)
pi_edna <- mb_copies / rowSums(mb_copies)
# Aliquots with zero copies across all target species produce NaN; fall
# back to uniform so rmultinom is well-defined (those rows yield few
# target reads anyway).
empty_rows <- !is.finite(rowSums(pi_edna))
pi_edna[empty_rows, ] <- 1 / n_species

# Multinomial draw of target reads per aliquot; stack as N_mb_long × n_species
mb_reads <- t(vapply(
  seq_len(N_mb_long),
  function(i) rmultinom(1, size = target_read_depth[i], prob = pi_edna[i, ])[, 1],
  integer(n_species)
))
storage.mode(mb_reads) <- "integer"
mb_total   <- as.integer(rowSums(mb_reads))    # = target_read_depth

# ---------------------------------------------------------------------------
# 10. Bundle and save
# ---------------------------------------------------------------------------
sim <- list(
  meta = list(
    n_species             = n_species,
    sp_names              = sp_names,
    sp_common             = sp_common,
    has_qpcr              = has_qpcr,
    n_qpcr_rep            = n_qpcr_rep,
    n_mb_rep              = n_mb_rep,
    prop_reduced_qpcr_rep = prop_reduced_qpcr_rep,
    prop_reduced_mb_rep   = prop_reduced_mb_rep,
    marker                = marker,
    n_stations            = n_stations,
    sample_depths         = sample_depths,
    n_sample_depth        = n_sample_depth,
    vol_filtered          = vol_filtered,
    vol_aliquot           = vol_aliquot,
    n_mb_reads            = n_mb_reads,
    N                     = N,
    N_qpcr_long           = N_qpcr_long,
    N_mb_long             = N_mb_long,
    utm_zone              = "10N",
    X_min_utm             = X_min,
    X_max_utm             = X_max,
    Y_min_utm             = Y_min,
    Y_max_utm             = Y_max,
    X_km_max              = X_km_max,
    Y_km_max              = Y_km_max,
    rotation_deg          = rotation_deg,
    bathy_profile         = list(
      rotation_deg       = rotation_deg,
      shelf_break_xprime = 180,
      abyssal_depth      = 2500,
      shelf_depth        =   80,
      slope_k            = 0.06
    ),
    coord_note = paste(
      "X (km): 0 = western edge of domain (100000 UTM E),",
      "500 = eastern edge (600000 UTM E); ",
      "Y (km): 0 = SF (~4180000 UTM N, 37.77°N),",
      "1270 = US/Canada border (~5450000 UTM N, ~49°N); ",
      "Z_bathy: bottom depth (m), derived from rotated X_prime; ",
      "Z_sample: water column sample depth (m)"
    )
  ),
  design = list(
    stations     = stations,
    samples      = samples,
    n_qpcr_rep_i = n_qpcr_rep_i,   # per-sample qPCR rep count (1, 2 or 3)
    n_mb_rep_i   = n_mb_rep_i      # per-sample MB rep count   (1, 2 or 3)
  ),
  truth = list(
    gp_field_si    = gp_field_si,    # f_s — includes habitat structure
    gp_mean_si     = gp_mean_si,     # diagnostic: GP mean used for sim
    lambda_true_si = lambda_true_si,
    C_obs_si       = C_obs_si,       # N × S  bottle copies
    qpcr_copies    = qpcr_copies,    # N_qpcr_long  aliquot copies (hake)
    mb_copies      = mb_copies,      # N_mb_long × S  aliquot copies
    zsample_effect = zsample_effect,
    pi_edna        = pi_edna,        # N_mb_long × S
    gp_params      = gp_params,
    qpcr_params    = qpcr_p
  ),
  observed = list(
    # qPCR long form (one row per qPCR aliquot)
    qpcr_sample_idx = qpcr_sample_idx,  # length N_qpcr_long
    qpcr_detect     = qpcr_detect,      # 0/1, length N_qpcr_long
    qpcr_ct         = qpcr_ct,          # NA when qpcr_detect == 0
    # Metabarcoding long form (one row per MB aliquot).
    # mb_reads / mb_total are target-species only; junk reads are stored
    # separately for transparency but are not passed to the Stan model.
    mb_sample_idx   = mb_sample_idx,    # length N_mb_long
    mb_reads        = mb_reads,         # N_mb_long × n_species
    mb_total        = mb_total,         # length N_mb_long (= target reads)
    mb_junk_reads   = junk_reads,       # length N_mb_long
    mb_pi_junk      = pi_junk           # length N_mb_long, 5–95%
  )
)

saveRDS(sim, "outputs/whale_edna_sim_v4.rds")

# ---------------------------------------------------------------------------
# 11. Detection-rate diagnostic
# ---------------------------------------------------------------------------
# A sample "has species s" under MB if any aliquot of that sample has
# reads[, s] > 0. A sample "has hake" under qPCR if any replicate of that
# sample has detect == 1. Target rates (real-data reference):
#   qPCR hake: 95%     MB hake: 95%
#   MB humpback: 20%   MB PWSD: 20%  (independent distributions)
cat("\n--- Detection rates (sample-level) ---\n")
qpcr_sample_det <- tapply(qpcr_detect, qpcr_sample_idx,
                          function(x) as.integer(any(x == 1)))
cat(sprintf("  qPCR  %-28s : %.1f%%\n", sp_common[1],
            100 * mean(qpcr_sample_det)))
for (s in seq_len(n_species)) {
  sample_det_s <- tapply(mb_reads[, s] > 0, mb_sample_idx, any)
  cat(sprintf("  MB    %-28s : %.1f%%\n",
              sp_common[s], 100 * mean(sample_det_s)))
}
cat("\n--- Mean / median lambda (animals/km^2) ---\n")
for (s in seq_len(n_species)) {
  cat(sprintf("  %-28s  mean=%7.3f  median=%7.3f\n",
              sp_common[s],
              mean(lambda_true_si[, s]),
              median(lambda_true_si[, s])))
}
