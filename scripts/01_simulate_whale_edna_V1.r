# =============================================================================
# simulate_whale_edna.R
#
# Anisotropic 3-species eDNA simulation
# Domain: Oregon / Washington coast, UTM Zone 10N
#
#   X (Easting)  : 200,000 – 500,000 m  (200 km, offshore to nearshore)
#   Y (Northing) : 4,800,000 – 5,200,000 m  (400 km, Cape Blanco to Strait)
#   Z_bathy      : bottom depth at each station (m), derived from X as a
#                  smooth function approximating the continental shelf/slope
#
# Two distinct depth processes:
#   1. Z_bathy — bathymetric depth at station (X, Y): drives TRUE animal
#      density via species habitat preferences (e.g. hake prefer 150–400 m
#      bottom depth; humpback prefer shallow shelf; PWSD prefer deep slope)
#      This is an intrinsic property of the location, not the sample.
#
#   2. Z_sample — depth at which water sample was collected: drives eDNA
#      vertical distribution within the water column. Animals may be at a
#      location but their eDNA is more concentrated at certain depths
#      (e.g. PWSD shed DNA near the surface; humpback during feeding dives).
#
# Species:
#   1. Merluccius productus  (Pacific hake)
#        — qPCR (3 replicates) + metabarcoding MARVER3 (3 replicates)
#        — shelf-slope break, prefers Z_bathy 150–400 m
#        — eDNA peaks at sample depths 150–300 m
#
#   2. Megaptera novaeangliae  (humpback whale)
#        — metabarcoding MARVER3 only (3 replicates)
#        — nearshore shelf, prefers Z_bathy < 200 m
#        — eDNA peaks at sample depth ~50–150 m (feeding dives)
#
#   3. Lagenorhynchus obliquidens  (Pacific white-sided dolphin)
#        — metabarcoding MARVER3 only (3 replicates)
#        — offshore slope, prefers Z_bathy > 500 m
#        — eDNA peaks at surface (0–50 m; active near-surface behaviour)
#
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
n_qpcr_rep <- 3
n_mb_rep   <- 3
marker     <- "MARVER3"

# UTM Zone 10N bounding box (metres)
X_min <- 200000;  X_max <- 500000   # 300 km E-W
Y_min <- 4800000; Y_max <- 5200000  # 400 km N-S

n_stations    <- 150

# Water column sample depths (Z_sample)
sample_depths <- c(0, 50, 150, 300, 500)   # m
n_sample_depth <- length(sample_depths)

vol_filtered  <- 2.0    # litres per sample
n_mb_reads    <- 5000   # MARVER3 reads per replicate

# ---------------------------------------------------------------------------
# 1. Station locations and bathymetry
#
# Bathymetric depth (Z_bathy) is approximated as a smooth function of
# Easting X, representing the continental shelf/slope profile off Oregon/WA:
#   - Nearshore (X > 250 km from offshore): shallow shelf, Z_bathy < 100 m
#   - Shelf break (~X = 180 km from offshore): Z_bathy ~ 200 m
#   - Slope (X = 80–180 km): Z_bathy 200–1500 m
#   - Offshore / abyssal (X < 80 km): Z_bathy > 1500 m
#
# We add moderate random noise to Z_bathy to reflect submarine canyons,
# banks, and real bathymetric variability.
# ---------------------------------------------------------------------------
set.seed(1)

# Stratified station design: denser on shelf/slope
n_shelf    <- round(n_stations * 0.65)
n_offshore <- n_stations - n_shelf

shelf_sta <- data.frame(
  X_utm = runif(n_shelf,    300000, 490000),
  Y_utm = runif(n_shelf,    Y_min,  Y_max)
)
offshore_sta <- data.frame(
  X_utm = runif(n_offshore, X_min,  320000),
  Y_utm = runif(n_offshore, Y_min,  Y_max)
)

stations <- bind_rows(shelf_sta, offshore_sta) %>%
  slice_sample(n = n_stations) %>%
  mutate(
    station = row_number(),
    X = (X_utm - X_min) / 1000,   # 0–300 km (0 = offshore, 300 = nearshore)
    Y = (Y_utm - Y_min) / 1000    # 0–400 km
  )

# Bathymetric depth as function of X (km from offshore boundary)
# Uses a logistic shelf-slope profile + random noise
bathy_mean_fn <- function(X_km) {
  # Shelf break at X ~ 180 km; abyssal depth ~ 2500 m; shelf ~ 80 m
  abyssal  <- 2500
  shelf    <- 80
  slope    <- abyssal + (shelf - abyssal) /
    (1 + exp(-0.06 * (X_km - 180)))
  pmax(slope, 10)   # minimum 10 m depth
}

set.seed(2)
stations <- stations %>%
  mutate(
    Z_bathy_mean = bathy_mean_fn(X),
    # Add spatially correlated noise (submarine canyons, banks)
    Z_bathy_noise = rnorm(n_stations, 0, Z_bathy_mean * 0.12),
    Z_bathy = pmax(Z_bathy_mean + Z_bathy_noise, 10)
  )

cat(sprintf("Stations: %d\n", n_stations))
cat(sprintf("  Z_bathy range: %.0f – %.0f m\n",
            min(stations$Z_bathy), max(stations$Z_bathy)))

# Full sample table: each station × each sample depth
# Only sample at depths <= Z_bathy (can't sample below the seafloor)
samples <- expand.grid(
  station     = 1:n_stations,
  depth_idx   = 1:n_sample_depth
) %>%
  left_join(stations, by = "station") %>%
  mutate(
    Z_sample  = sample_depths[depth_idx],
    sample_id = row_number()
  ) %>%
  # Drop samples where Z_sample > Z_bathy (below seafloor)
  filter(Z_sample <= Z_bathy) %>%
  mutate(sample_id = row_number())   # reindex after filtering

N <- nrow(samples)
cat(sprintf("Samples after filtering Z_sample > Z_bathy: %d\n", N))

# GP coordinates use (X, Y, Z_bathy) — bathymetry drives true density
coords_gp <- as.matrix(samples[, c("X", "Y", "Z_bathy")])

# ---------------------------------------------------------------------------
# 2. GP hyperparameters
#
# The GP field is a function of (X, Y, Z_bathy).
# Z_bathy length-scale reflects how sharply density changes across the
# shelf/slope bathymetric gradient.
#
# Hake:      shelf-slope break specialist; Z_bathy lz ~ 150 m
#            (density changes over ~150 m of bottom depth)
# Humpback:  shallow shelf specialist; Z_bathy lz ~ 100 m
#            (sharp preference for shallow banks and shelf)
# PWSD:      deep slope / oceanic; Z_bathy lz ~ 300 m
#            (broadly distributed across slope depths)
#
# Z_bathy mean function: Gaussian-shaped preference over Z_bathy values
# rather than just X, making the depth preference explicit.
# ---------------------------------------------------------------------------

gp_params <- list(
  
  hake = list(
    name  = sp_names[1],
    sigma = 1.2,
    lx    =  50,    # km — cross-shore
    ly    = 150,    # km — along-coast
    lz    = 150,    # m  — bathymetric gradient
    mu    = log(20),
    # Bathymetric preference: peak at Z_bathy ~ 200 m (shelf-slope break)
    zbathy_pref_mu = 200,   # m
    zbathy_pref_sd = 100,   # m
    zbathy_pref_amp = 2.0,  # log-scale amplitude of preference
    # eDNA vertical distribution in water column (Z_sample offsets)
    # Hake live and shed DNA at mid-depth; eDNA peaks at 150–300 m
    zsample_pref = c(-2.0, -0.8,  1.0,  0.8, -0.5)  # at 0,50,150,300,500 m
  ),
  
  humpback = list(
    name  = sp_names[2],
    sigma = 1.0,
    lx    =  40,
    ly    = 200,
    lz    = 100,
    mu    = log(5),
    # Bathymetric preference: nearshore shelf, Z_bathy 50–150 m
    zbathy_pref_mu = 100,
    zbathy_pref_sd =  80,
    zbathy_pref_amp = 2.0,
    # eDNA water column: peaks at 50–150 m (feeding dive depths)
    # Detectable at surface; drops off at 300–500 m
    zsample_pref = c(-0.5,  0.8,  0.6, -0.5, -1.8)
  ),
  
  pwsd = list(
    name  = sp_names[3],
    sigma = 1.3,
    lx    =  30,
    ly    = 100,
    lz    = 300,
    mu    = log(8),
    # Bathymetric preference: deep slope / oceanic, Z_bathy > 500 m
    zbathy_pref_mu = 800,
    zbathy_pref_sd = 400,
    zbathy_pref_amp = 2.0,
    # eDNA water column: surface-active; eDNA concentrated near surface
    zsample_pref = c( 1.5,  0.8, -0.5, -1.0, -2.0)
  )
)

# ---------------------------------------------------------------------------
# 3. Draw GP realisations
#
# log(lambda_si) = mu_s
#                + zbathy_pref(Z_bathy_i)    # bathymetric habitat effect
#                + f_s(X_i, Y_i, Z_bathy_i) # spatially structured residual
#
# Note: Z_bathy enters the GP kernel, not Z_sample.
# The zbathy_pref mean function is separate from the GP kernel and encodes
# the species-specific bathymetric habitat association directly.
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

zbathy_pref_fn <- function(Z_bathy, mu_z, sd_z, amp) {
  # Gaussian-shaped bathymetric preference, scaled to [-amp/2, +amp/2]
  raw <- dnorm(Z_bathy, mu_z, sd_z)
  raw <- raw / max(raw)        # 0 to 1
  amp * (raw - 0.5)            # centred: peak = +amp/2, trough = -amp/2
}

cat("Drawing GP realisations (N =", N, ")...\n")
lambda_true      <- matrix(NA, N, n_species)
gp_field         <- matrix(NA, N, n_species)
zbathy_mean_field <- matrix(NA, N, n_species)

for (s in seq_len(n_species)) {
  p  <- gp_params[[s]]
  
  # GP kernel over (X, Y, Z_bathy)
  K  <- aniso_cov(coords_gp, p$sigma, p$lx, p$ly, p$lz)
  fs <- as.vector(mvrnorm(1, rep(0, N), K))
  gp_field[, s] <- fs
  
  # Bathymetric habitat preference (mean function)
  z_mu <- zbathy_pref_fn(samples$Z_bathy,
                         p$zbathy_pref_mu,
                         p$zbathy_pref_sd,
                         p$zbathy_pref_amp)
  zbathy_mean_field[, s] <- z_mu
  
  # True animal density: function of X, Y, Z_bathy only
  # (before water column sampling process)
  lambda_true[, s] <- exp(p$mu + z_mu + fs)
}

# ---------------------------------------------------------------------------
# 4. eDNA water column process
#
# The observed eDNA concentration at sample depth Z_sample is:
#   C_obs_si ~ ZI-Poisson(lambda_si * zsample_effect_si * vol)
#
# zsample_effect: multiplicative factor on eDNA concentration based on
# where in the water column the sample was taken relative to where the
# animal sheds DNA. This is interpolated from the zsample_pref offsets.
#
# This separates the ecological question (where are the animals, as a
# function of X/Y/Z_bathy) from the sampling question (where in the water
# column is eDNA detectable).
# ---------------------------------------------------------------------------

# Interpolated zsample effect for each sample
zsample_effect <- matrix(NA, N, n_species)

for (s in seq_len(n_species)) {
  p <- gp_params[[s]]
  # Interpolate zsample_pref over the actual sample_depths
  pref_fn <- approxfun(
    x      = sample_depths,
    y      = p$zsample_pref,
    method = "linear",
    rule   = 2
  )
  zsample_effect[, s] <- exp(pref_fn(samples$Z_sample))
}

# Effective eDNA concentration = lambda (animal density) * water column factor
lambda_edna <- lambda_true * zsample_effect   # N × S

# ---------------------------------------------------------------------------
# 5. Zero-inflated Poisson: copies in filtered volume
# ---------------------------------------------------------------------------
p_zi <- c(0.30, 0.45, 0.50)

C_obs   <- matrix(NA_integer_, N, n_species)
zi_flag <- matrix(NA_integer_, N, n_species)

for (s in seq_len(n_species)) {
#  zi           <- rbinom(N, 1, p_zi[s]) # is zi the bottle or the aliquot?
  # if it's the aliquot, there's no zero inflation
  # if it's the bottle, there could be zero-inflation, but it should be a function
  # of the underlying density, not some arbitrary number
  pois_draw    <- rpois(N, lambda_edna[, s] * vol_filtered)
  C_obs[, s]   <- pois_draw
  zi_flag[, s] <- ifelse(pois_draw == 0, 0, 1)
}
# note that qPCRs are drawn from the bottle, metabarcoding is separate
# need to add a layer for aliquots within bottles
# ---------------------------------------------------------------------------
# 6. qPCR hurdle — hake only (3 replicates)
# ---------------------------------------------------------------------------
qpcr_p <- list(
  kappa    = 0.85,
  alpha_ct = 38.0,
  beta_ct  = 3.32,
  sigma_ct = 0.50
)

qpcr_detect_raw <- matrix(NA_integer_, N, n_qpcr_rep)
qpcr_ct_raw     <- matrix(NA_real_,    N, n_qpcr_rep)

p_det_hake <- 1 - exp(-qpcr_p$kappa * C_obs[, 1])

for (r in seq_len(n_qpcr_rep)) {
  det <- rbinom(N, 1, p_det_hake)
  Ct  <- ifelse(
    det == 1L,
    qpcr_p$alpha_ct -
      qpcr_p$beta_ct * log10(pmax(C_obs[, 1], 1)) +
      rnorm(N, 0, qpcr_p$sigma_ct),
    NA_real_
  )
  qpcr_detect_raw[, r] <- det
  qpcr_ct_raw[, r]     <- Ct
}

qpcr_n_detect   <- rowSums(qpcr_detect_raw)
qpcr_any_detect <- as.integer(qpcr_n_detect > 0)
qpcr_mean_ct    <- apply(qpcr_ct_raw, 1,
                         function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))

# ---------------------------------------------------------------------------
# 7. Metabarcoding Beta-Binomial — all species, 3 MARVER3 replicates
# ---------------------------------------------------------------------------
phi_bb  <- c(10.0, 7.0, 7.0)
# IMPORTANT
# phi_bb should be a function of the concentration
# there is a formula for this
# this should be a multinomial 

# Proportions based on eDNA concentration (lambda_edna, not lambda_true)
# because metabarcoding reflects what is in the water at sampling depth

# this needs to be linked to the same bottle as the qPCR replicates, 
# not directly from the density surface
pi_edna <- lambda_edna / rowSums(lambda_edna)

rbetabinom <- function(n, size, mu, phi) {
  mu  <- pmin(pmax(mu,  1e-6), 1 - 1e-6)
  phi <- pmax(phi, 0.1)
  p   <- rbeta(n, mu * phi, (1 - mu) * phi)
  rbinom(n, size, p)
}

mb_reads_rep <- array(NA_integer_, dim = c(N, n_species, n_mb_rep))
for (s in seq_len(n_species)) {
  for (rep in seq_len(n_mb_rep)) {
    mb_reads_rep[, s, rep] <- rbetabinom(
      N, size = n_mb_reads, # in reality n_mb_reads is stochastic, from data
      # simulation should draw a random number of total reads 10K-100K
      # total target reads might be 1000-10K
      # Pedro to check actual numbers, but mm usually around 1000
      # 1.4% of reads are MM
      mu = pi_edna[, s], phi = phi_bb[s]
    )
  }
}

mb_total_rep  <- apply(mb_reads_rep, c(1, 3), sum)
mb_reads_mean <- apply(mb_reads_rep, c(1, 2), mean)
mb_total_mean <- rowSums(mb_reads_mean)

# ---------------------------------------------------------------------------
# 8. Bundle and save
# ---------------------------------------------------------------------------
sim <- list(
  meta = list(
    n_species      = n_species,
    sp_names       = sp_names,
    sp_common      = sp_common,
    has_qpcr       = has_qpcr,
    n_qpcr_rep     = n_qpcr_rep,
    n_mb_rep       = n_mb_rep,
    marker         = marker,
    n_stations     = n_stations,
    sample_depths  = sample_depths,
    n_sample_depth = n_sample_depth,
    vol_filtered   = vol_filtered,
    n_mb_reads     = n_mb_reads,
    N              = N,
    utm_zone       = "10N",
    X_min_utm      = X_min,
    X_max_utm      = X_max,
    Y_min_utm      = Y_min,
    Y_max_utm      = Y_max,
    coord_note     = paste(
      "X: 0=offshore (200000 m Easting), 300=nearshore (500000 m Easting);",
      "Y: 0=south (Cape Blanco ~4800000 m Northing), 400=north;",
      "Z_bathy: bottom depth at station (m);",
      "Z_sample: water column sample depth (m)"
    )
  ),
  design = list(
    stations = stations,   # one row per station: X, Y, Z_bathy, X_utm, Y_utm
    samples  = samples     # one row per sample: station, X, Y, Z_bathy,
    #   Z_sample, depth_idx, sample_id
  ),
  truth = list(
    gp_field          = gp_field,          # N × S  (GP residual only)
    zbathy_mean_field = zbathy_mean_field,  # N × S  (bathymetric mean fn)
    lambda_true       = lambda_true,        # N × S  (true animal density)
    lambda_edna       = lambda_edna,        # N × S  (eDNA conc at Z_sample)
    zsample_effect    = zsample_effect,     # N × S  (water column multiplier)
    pi_edna           = pi_edna,            # N × S  (eDNA proportions)
    zi_flag           = zi_flag,            # N × S
    gp_params         = gp_params,
    qpcr_params       = qpcr_p,
    p_zi              = p_zi,
    phi_bb            = phi_bb
  ),
  observed = list(
    # qPCR (hake only)
    qpcr_detect_raw = qpcr_detect_raw,   # N × n_qpcr_rep
    qpcr_ct_raw     = qpcr_ct_raw,       # N × n_qpcr_rep
    qpcr_n_detect   = qpcr_n_detect,     # N
    qpcr_any_detect = qpcr_any_detect,   # N
    qpcr_mean_ct    = qpcr_mean_ct,      # N  (NA if no detection)
    # Metabarcoding (all species, MARVER3)
    mb_reads_rep    = mb_reads_rep,      # N × S × n_mb_rep
    mb_total_rep    = mb_total_rep,      # N × n_mb_rep
    mb_reads_mean   = mb_reads_mean,     # N × S  (diagnostic)
    mb_total_mean   = mb_total_mean      # N      (diagnostic)
  )
)

saveRDS(sim, "outputs/whale_edna_sim_V1.rds")

# ---------------------------------------------------------------------------
# 9. Summary printout
# ---------------------------------------------------------------------------
cat("\nSaved outputs/whale_edna_sim_V1.rds\n")
cat(sprintf("  Domain    : UTM Zone 10N\n"))
cat(sprintf("              X = [%.0f, %.0f] m Easting\n",  X_min, X_max))
cat(sprintf("              Y = [%.0f, %.0f] m Northing\n", Y_min, Y_max))
cat(sprintf("  Stations  : %d\n", n_stations))
cat(sprintf("  Z_bathy   : %.0f – %.0f m (mean %.0f m)\n",
            min(stations$Z_bathy), max(stations$Z_bathy),
            mean(stations$Z_bathy)))
cat(sprintf("  Samples   : %d (after removing Z_sample > Z_bathy)\n", N))
cat(sprintf("  Marker    : %s (%d replicates/sample)\n", marker, n_mb_rep))
cat("\n  Per-species summary:\n")

for (s in seq_len(n_species)) {
  cat(sprintf("\n  [%d] %s (%s)\n", s, sp_common[s], sp_names[s]))
  cat(sprintf("      True lambda  : mean = %.1f, range = %.1f – %.1f copies/L\n",
              mean(lambda_true[, s]),
              min(lambda_true[, s]),
              max(lambda_true[, s])))
  cat(sprintf("      eDNA lambda  : mean = %.1f (after water column effect)\n",
              mean(lambda_edna[, s])))
  cat(sprintf("      ZI           : %.0f%%\n", 100 * p_zi[s]))
  cat(sprintf("      Z_bathy pref : peak at %.0f m (SD %.0f m)\n",
              gp_params[[s]]$zbathy_pref_mu,
              gp_params[[s]]$zbathy_pref_sd))
  cat(sprintf("      Z_sample pref: %s\n",
              paste(sprintf("%.1f", gp_params[[s]]$zsample_pref),
                    collapse = " / ")))
  if (has_qpcr[s]) {
    cat(sprintf("      qPCR detect  : %.0f%% of samples\n",
                100 * mean(qpcr_any_detect)))
  } else {
    cat("      qPCR         : not collected\n")
  }
  mb_det <- mean(mb_reads_mean[, s] > 0)
  cat(sprintf("      MB detection : %.0f%% of samples\n", 100 * mb_det))
}

