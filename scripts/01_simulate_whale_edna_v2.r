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
#        — qPCR (3 replicates) + metabarcoding MARVER1 (3 replicates)
#        — shelf-slope break, prefers Z_bathy 150–400 m
#        — eDNA peaks at sample depths 150–300 m
#
#   2. Megaptera novaeangliae  (humpback whale)
#        — metabarcoding MARVER1 only (3 replicates)
#        — nearshore shelf, prefers Z_bathy < 200 m
#        — eDNA peaks at sample depth ~50–150 m (feeding dives)
#
#   3. Lagenorhynchus obliquidens  (Pacific white-sided dolphin)
#        — metabarcoding MARVER1 only (3 replicates)
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
R <- n_qpcr_rep # might want to use the specific things
marker     <- "MARVER1"

# UTM Zone 10N bounding box (metres)
X_min <- 200000;  X_max <- 500000   # 300 km E-W
Y_min <- 4800000; Y_max <- 5200000  # 400 km N-S

n_stations    <- 150
conv_factor <- 10 # conversion between animals and copies
# Water column sample depths (Z_sample)
sample_depths <- c(0, 50, 150, 300, 500)   # m
n_sample_depth <- length(sample_depths)

vol_filtered  <- 2.5    # litres per sample
vol_aliquot <- 2 # microliters per aliquot
n_mb_reads    <- 50000   # MARVER1 reads per replicate 

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
# could plug in actual bathymetry here
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
    zsample_pref = c(-2.0+2, -0.8+2,  1.0+2,  0.8+2, -0.5+2)  # at 0,50,150,300,500 m
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
    zsample_pref = c(-0.5+0.5,  0.8+0.5,  0.6+0.5, -0.5+0.5, -1.8+0.5)
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
    zsample_pref = c( 1.5-1.5,  0.8-1.5, -0.5-1.5, -1.0-1.5, -2.0-1.5)
  )
)

# ---------------------------------------------------------------------------
# 3. Draw GP realisations
# lambda_true is in units of animals/km^2
# log(lambda_true_si) = # lambda_true_si is species density at surface loc i
#               mu_s 
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
lambda_true_si      <- matrix(NA, N, n_species) # note doesn't match index order
gp_field_si         <- matrix(NA, N, n_species)
zbathy_mean_field_si <- matrix(NA, N, n_species)

for (s in seq_len(n_species)) {
  p  <- gp_params[[s]]
  
  # GP kernel over (X, Y, Z_bathy)
  K  <- aniso_cov(coords_gp, p$sigma, p$lx, p$ly, p$lz)
  fs <- as.vector(mvrnorm(1, rep(0, N), K))
  gp_field_si[, s] <- fs
  
  # Bathymetric habitat preference (mean function)
  z_mu <- zbathy_pref_fn(samples$Z_bathy,
                         p$zbathy_pref_mu,
                         p$zbathy_pref_sd,
                         p$zbathy_pref_amp)
  zbathy_mean_field_si[, s] <- z_mu
  
  # True animal density: function of X, Y, Z_bathy only
  # (before water column sampling process)
  lambda_true_si[, s] <- exp(p$mu + fs) # note could add z_mu back in here
} # end loop over species

# ---------------------------------------------------------------------------
# 4. eDNA water column process
#
# The observed (bottle) eDNA concentration at sample depth Z_sample is:
# For now, say 1:10 relationship between animals and concentration
# conv_factor = e.g., 10 copies/animal*L*km^2
#   C_obs_si ~ NB(conv_factor * lambda_true_si * zsample_effect * vol_filtered, 
#             phi_nb = 10)
# C_obs_si is in units of copies in the bottle
# Need to divide by volume eluted 100uL
# concentration in the aliquot is
# C_obs_si_a ~ Bin(C_obs_si, vol_aliquot/100uL)

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

#   C_obs_si ~ NB(conv_factor * lambda_true_si * zsample_effect_si * vol_filtered, 
#             theta = 10)

C_obs_si <- matrix(NA, N, n_species)
C_obs_si_a_mb <- matrix(NA, N*R, n_species+1)
C_obs_si_a_qpcr <- matrix(NA, N*R, n_species+1)

C_obs_si_a_mb[, n_species+1] <- rep(1:N, each = 3) # track which sample 
C_obs_si_a_qpcr[, n_species+1] <- rep(1:N, each = 3)

for (s in 1:n_species){
  
  # copies in the bottle
  C_obs_si[, s] <- MASS::rnegbin(N, mu = conv_factor * lambda_true_si[, s] * zsample_effect[, s] * vol_filtered,
                          theta = 10)
  # copies in the aliquot
  C_obs_si_a_mb[ ,s] <- rbinom(N*R, rep(C_obs_si[,s], each = R), vol_aliquot/100)
  if(s==1){C_obs_si_a_qpcr[,s] <- rbinom(N*R, rep(C_obs_si[,s], each = R), vol_aliquot/100)
  }} # end for s in n_species

# ---------------------------------------------------------------------------
# 6. qPCR hurdle — hake only (3 replicates)
# ---------------------------------------------------------------------------
qpcr_p <- list( 
  kappa    = 0.85,
  alpha_ct = 38.0,
  beta_ct  = 1.44, # outputs (conc) will be in log e
  sigma_ct = 0.50 # this assumes no variation with concentration, might want to change
)

qpcr_detect_raw <- matrix(NA_integer_, N*R, n_species)
qpcr_ct_raw     <- matrix(NA_real_,    N*R, n_species)

p_det_hake <- 1 - exp(-qpcr_p$kappa * C_obs_si_a_qpcr[, 1]) # hake only right now

  det <- rbinom(N*R, 1, p_det_hake)
  Ct  <- ifelse(
    det == 1L,
    qpcr_p$alpha_ct -
      qpcr_p$beta_ct * log(pmax(C_obs_si_a_qpcr[, 1], 1)) +
      rnorm(N*R, 0, qpcr_p$sigma_ct),
    NA_real_
  )
  qpcr_detect_raw[, 1] <- det
  qpcr_ct_raw[, 1]     <- Ct

# ---------------------------------------------------------------------------
# 7. Metabarcoding Beta-Binomial — all species, 3 MARVER1 replicates
# ---------------------------------------------------------------------------
# Per-aliquot multinomial draw: for row i, draw read_depth[i] reads from
# the species-proportion vector pi_edna[i, ]. rmultinom is not vectorised
# over a varying `size` with a vector prob, so we loop and stack.
# Result shape: n_species × (N*R), i.e. one column per sample-aliquot.

read_depth <- round(runif(N*R, n_mb_reads/2, n_mb_reads))

pi_edna <- C_obs_si_a_mb[, 1:n_species] / rowSums(C_obs_si_a_mb[, 1:n_species])

mb_reads_rep <- vapply(
  seq_len(N * R),
  function(i) rmultinom(1, size = read_depth[i], prob = pi_edna[i, ])[, 1],
  integer(n_species)
)


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
    gp_field_si          = gp_field_si,          # N × S  (GP residual only)
    zbathy_mean_field_si = zbathy_mean_field_si,  # N × S  (bathymetric mean fn)
    lambda_true_si       = lambda_true_si,        # N × S  (true animal density)
    C_obs_si = C_obs_si,
    C_obs_si_a_mb = C_obs_si_a_mb,
    C_obs_si_a_qpcr = C_obs_si_a_mb,
    zsample_effect    = zsample_effect,     # N × S  (water column multiplier)
    pi_edna           = pi_edna,            # N × S  (eDNA proportions)
    gp_params         = gp_params,
    qpcr_params       = qpcr_p
    ),
  observed = list(
    # qPCR (hake only)
    qpcr_detect_raw = qpcr_detect_raw,   # N × n_qpcr_rep
    qpcr_ct_raw     = qpcr_ct_raw,       # N × n_qpcr_rep
    # Metabarcoding (all species, MARVER1)
    mb_reads_rep    = mb_reads_rep      # N × S × n_mb_rep
  )
)

saveRDS(sim, "outputs/whale_edna_sim_v2.rds")
