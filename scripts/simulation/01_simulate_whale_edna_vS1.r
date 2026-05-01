# =============================================================================
# simulate_whale_edna_vS1.R
#
# vS1 — simplified 2-D simulation. Whale density is a function of (X, Y)
# only; Z_bathy and the rotated cross-shore axis are removed entirely.
# Stations are drawn uniformly from the study area (no shelf / slope /
# offshore stratification). All sampling is surface-only.
#
# Latent field: zero-mean 2-D anisotropic GP
#     log(lambda_si) = mu_s + f_s(X, Y),    f_s ~ GP(0, K(lx, ly))
#
# qPCR calibration parameters (kappa, sigma_ct, alpha_ct, beta_ct) are
# all treated as known in downstream stages (sigma_ct passed inflated
# to absorb discrete-count noise; see scripts/03_format_stan_data_v4.2.r).
#
# Includes changes from v4.1:
#   * collaborator changes to the Ct sigma function
#   * collaborator reformulation of the junk category so junk reads are
#     retained in mb_reads
# =============================================================================

library(tidyverse)
library(MASS)
set.seed(20260501)

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
prop_reduced_qpcr_rep <- 0 # set to 0 for this simulation case
prop_reduced_mb_rep   <- 0

# UTM Zone 10N bounding box (metres): San Francisco (~37.77°N) to 49°N
X_min <- 100000;  X_max <- 600000   # 500 km E-W
Y_min <- 4180000; Y_max <- 5450000  # ~1270 km N-S

# Domain extent in km (for all GP / preference calculations)
X_km_max <- (X_max - X_min) / 1000   # 500 km
Y_km_max <- (Y_max - Y_min) / 1000   # 1270 km

n_stations    <- 200

# Species-specific animals -> eDNA-copies conversion factor (copies/L per
# animal/km^2 per litre filtered). Whales shed far more eDNA per animal
# than hake, so conv_factor is much larger for the two whale species.
# Ordered to match gp_params: hake, humpback, PWSD. Same values as v4.
conv_factor <- c(hake = 10, humpback = 200, pwsd = 110, junk=10)

# Surface-only sampling: a single sample depth at Z_sample = 0.
sample_depths <- c(0)
n_sample_depth <- length(sample_depths)

vol_filtered  <- 2.5
vol_aliquot <- 2

# MB read-depth distribution per aliquot (matches real sequencing runs).
# A two-component log-normal mixture: 80% "typical" reads centred on 75k,
# with 20% from a wider "problem-run" component that can dip to a few
# thousand or spike past 200k. Hard-clipped to the observed real-data
# range of [1000, 250000].
mb_reads_p_tight       <- 0.80
mb_reads_tight_meanlog <- log(75000)
mb_reads_tight_sdlog   <- 0.25
mb_reads_wide_meanlog  <- log(40000)
mb_reads_wide_sdlog    <- 1.00
mb_reads_min           <- 1000L
mb_reads_max           <- 250000L

# ---------------------------------------------------------------------------
# 1. Station locations
#
# Uniform random sampling over the (X, Y) study area — no rotation, no
# bathymetry-based stratification. Z_bathy is removed entirely; the GP
# below is 2-D over (X, Y).
# ---------------------------------------------------------------------------
set.seed(1)

stations <- tibble(
  station = seq_len(n_stations),
  X       = runif(n_stations, 0, X_km_max),   # km E-W from western edge
  Y       = runif(n_stations, 0, Y_km_max),   # km N-S from SF
  X_utm   = X_min + X * 1000,
  Y_utm   = Y_min + Y * 1000
)

cat(sprintf("Stations: %d (uniform over %.0f x %.0f km domain)\n",
            n_stations, X_km_max, Y_km_max))

# Full sample table: station × sample depth. Surface-only here, so this
# is just one row per station; the depth_idx loop is preserved for
# shape compatibility with downstream stages that expect a `samples`
# table indexed by both station and depth_idx.
samples <- expand.grid(
  station   = 1:n_stations,
  depth_idx = 1:n_sample_depth
) %>%
  left_join(stations, by = "station") %>%
  mutate(
    Z_sample  = sample_depths[depth_idx],
    sample_id = row_number()
  )

N <- nrow(samples)
cat(sprintf("Samples (surface-only): %d\n", N))

# GP coordinates: (X, Y) only. The latent field is 2-D.
coords_gp <- as.matrix(samples[, c("X", "Y")])

# ---------------------------------------------------------------------------
# 2. GP hyperparameters
#
# Model form exposed to downstream inference:
#   log(lambda_s) = mu_s + f_s(X, Y)
#
# 2-D zero-mean anisotropic GP per species, with cross-shore (lx) and
# along-shore (ly) length-scales. Z_bathy is no longer part of the field.
# `zsample_pref` is kept for shape compatibility with downstream code,
# but with surface-only sampling it has a single entry (no depth effect).
# ---------------------------------------------------------------------------

gp_params <- list(

  hake = list(
    name  = sp_names[1],
    sigma = 2,
    lx    =  50,
    ly    = 300,
    mu    = log(15),
    # Surface-only sampling: zsample_pref aligns with sample_depths = c(0)
    zsample_pref = c(0.0)
  ),

  humpback = list(
    name  = sp_names[2],
    sigma = 1.0,
    lx    =  50,
    ly    = 300,
    mu    = log(2.5), # made larger than is realistic
    zsample_pref = c(0.0)
  ),

  pwsd = list(
    name  = sp_names[3],
    sigma = 1.3,
    lx    =  40,
    ly    = 300,
    mu    = log(5), # made larger than is realistic
    zsample_pref = c(0.0)
  )
)

# ---------------------------------------------------------------------------
# 3. GP draws (2-D)
# ---------------------------------------------------------------------------

aniso_cov <- function(coords, sigma, lx, ly) {
  n <- nrow(coords)
  K <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in i:n) {
      d2 <- ((coords[i,1]-coords[j,1])/lx)^2 +
            ((coords[i,2]-coords[j,2])/ly)^2
      K[i,j] <- K[j,i] <- sigma^2 * exp(-0.5 * d2)
    }
  }
  K + diag(1e-6, n)
}

cat("Drawing zero-mean GP realisations (N =", N, ")...\n")
lambda_true_si <- matrix(NA, N, n_species)
gp_field_si    <- matrix(NA, N, n_species)
gp_mean_si     <- matrix(0,  N, n_species)   # zero-mean - kept for output shape

for (s in seq_len(n_species)) {
  p  <- gp_params[[s]]

  K  <- aniso_cov(coords_gp, p$sigma, p$lx, p$ly)

  # Zero-mean GP: f_s ~ MVN(0, K_s).
  fs <- as.vector(mvrnorm(1, mu = rep(0.0, N), Sigma = K))
  gp_field_si[, s] <- fs

  lambda_true_si[, s] <- exp(p$mu + fs)
}

# ---------------------------------------------------------------------------
# 5. eDNA water column process
# ---------------------------------------------------------------------------

zsample_effect <- matrix(NA, N, n_species)
for (s in seq_len(n_species)) {
  p <- gp_params[[s]]
  if (n_sample_depth == 1L) {
    # Surface-only: every sample takes the single zsample_pref value.
    zsample_effect[, s] <- exp(p$zsample_pref[1])
  } else {
    pref_fn <- approxfun(
      x      = sample_depths,
      y      = p$zsample_pref,
      method = "linear",
      rule   = 2
    )
    zsample_effect[, s] <- exp(pref_fn(samples$Z_sample))
  }
}

C_obs_si <- matrix(NA_integer_, N, n_species)
for (s in 1:n_species) {
  C_obs_si[, s] <- MASS::rnegbin(
    N,
    mu    = conv_factor[s] * lambda_true_si[, s] * zsample_effect[, s] * vol_filtered,
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

qpcr_exp_copies <- C_obs_si[qpcr_sample_idx, 1] * (vol_aliquot / 100)

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
  gamma_0  = 0.30,  
  gamma_1  = -0.30,
  sigma_0 = 0.23
)

p_det_hake   <- 1 - exp(-qpcr_p$kappa * qpcr_exp_copies)
qpcr_detect  <- rbinom(N_qpcr_long, 1, p_det_hake)
qpcr_p$sigma_ct     <- sqrt(qpcr_p$sigma_0^2 + 
                       exp(2*(qpcr_p$gamma_0 + qpcr_p$gamma_1 * log(qpcr_exp_copies))))
qpcr_ct      <- ifelse(
  qpcr_detect == 1L,
  qpcr_p$alpha_ct -
    qpcr_p$beta_ct * log(pmax(qpcr_exp_copies, 1)) +
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

# For every metabarcoding sample, define an expected copies of junk
# and then a realized number of copies in each replicate.
# Junk does not have a spatial pattern.
junk_mean <- 50
junk_theta  <- 1000
junk_exp_copies <- rnegbin(nrow(samples),mu=junk_mean,theta=junk_theta)

C_obs_si_junk <- junk_exp_copies * 50

mb_junk_copies <- matrix(NA_integer_, N_mb_long, 1)
mb_junk_copies <- rbinom(N_mb_long,
                           size = C_obs_si_junk[mb_sample_idx],
                           prob = vol_aliquot / 100)

# Per-aliquot total read depth from the mixture defined at the top of
# the script.
# Per-aliquot total read depth from the mixture defined at the top of
# the script.
component_tight <- rbinom(N_mb_long, 1, mb_reads_p_tight)
read_depth <- ifelse(
  component_tight == 1L,
  rlnorm(N_mb_long, mb_reads_tight_meanlog, mb_reads_tight_sdlog),
  rlnorm(N_mb_long, mb_reads_wide_meanlog,  mb_reads_wide_sdlog)
)
read_depth <- pmin(pmax(as.integer(round(read_depth)), mb_reads_min), mb_reads_max)

pi_junk    <- mb_junk_copies / (rowSums(mb_copies) + mb_junk_copies)
# hist(pi_junk)

# Target-species proportions (within the non-junk pool)
pi_edna <- mb_copies / (rowSums(mb_copies) + mb_junk_copies)
pi_all <- cbind(pi_edna,pi_junk)

# Aliquots with zero copies across all target species have no target
# signal at all — all reads in that aliquot are junk. Force
# target_read_depth = 0 for those rows so the multinomial produces zero
# target reads, and bump the junk read count to preserve the per-aliquot
# total. A sentinel pi is set only to keep rmultinom well-defined when
# size = 0.

empty_rows <- !is.finite(rowSums(pi_edna))
pi_edna[empty_rows, ] <- 0
#target_read_depth[empty_rows] 

# Multinomial draw of target reads per aliquot; stack as N_mb_long × n_species
mb_reads <- t(vapply(
  seq_len(N_mb_long),
  function(i) rmultinom(1, size = read_depth[i], prob = pi_all[i, ])[, 1],
  integer(n_species+1)
))
storage.mode(mb_reads) <- "integer"
mb_total   <- as.integer(rowSums(mb_reads))    # = target_read_depth

target_read_depth <- read_depth - mb_reads[,"pi_junk"]

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
    conv_factor           = conv_factor,     # named vector, length n_species
    mb_reads_distribution = list(
      p_tight       = mb_reads_p_tight,
      tight_meanlog = mb_reads_tight_meanlog,
      tight_sdlog   = mb_reads_tight_sdlog,
      wide_meanlog  = mb_reads_wide_meanlog,
      wide_sdlog    = mb_reads_wide_sdlog,
      min           = mb_reads_min,
      max           = mb_reads_max
    ),
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
    coord_note = paste(
      "X (km): 0 = western edge of domain (100000 UTM E),",
      "500 = eastern edge (600000 UTM E); ",
      "Y (km): 0 = SF (~4180000 UTM N, 37.77°N),",
      "1270 = US/Canada border (~5450000 UTM N, ~49°N); ",
      "Z_sample: water column sample depth (m, surface only)"
    )
  ),
  design = list(
    stations     = stations,
    samples      = samples,
    n_qpcr_rep_i = n_qpcr_rep_i,   # per-sample qPCR rep count (1, 2 or 3)
    n_mb_rep_i   = n_mb_rep_i      # per-sample MB rep count   (1, 2 or 3)
  ),
  truth = list(
    gp_field_si    = gp_field_si,    # f_s — zero-mean 2-D GP draw over (X, Y)
    gp_mean_si     = gp_mean_si,     # zero matrix; kept for output shape
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
    mb_reads        = mb_reads,         # N_mb_long × n_species + 1
    mb_total        = mb_total,         # length N_mb_long (= target reads)
    mb_junk_reads   = mb_reads[,ncol(mb_reads)],       # length N_mb_long
    mb_pi_junk      = pi_junk           # length N_mb_long, 5–95%
  )
)

dir.create("outputs/whale_edna_output_v4.2", showWarnings = FALSE, recursive = TRUE)
saveRDS(sim, "outputs/whale_edna_output_v4.2/whale_edna_sim_v4.2.rds")

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
