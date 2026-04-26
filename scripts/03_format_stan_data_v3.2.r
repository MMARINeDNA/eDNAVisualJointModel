# =============================================================================
# 03_format_stan_data_v3.2.R
#
# Read the simulated eDNA data (outputs/whale_edna_output_v3.2/whale_edna_sim_v3.2.rds)
# and assemble the Stan data list used by stan/whale_edna_hsgp_v3.2.stan. Writes:
#
#   outputs/whale_edna_output_v3.2/stan_data.rds
#
# Differences from v3:
#   * S = 1 (single species, hake). Surface-only samples.
#   * Metabarcoding inputs removed entirely (the v3.2 Stan model has no
#     MB likelihood).
#   * `kappa` is passed as data (fixed), not sampled.
#   * MB-only prior hyperparameters dropped.
# =============================================================================

library(tidyverse)

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------

# HSGP basis terms per dimension. Riutort-Mayol et al. (2023) rule of
# thumb: M >= 3.2 * c * (S / l), where c = boundary factor (1.5),
# S = data half-range in raw units (250 km, 635 km, 1750 m after the
# v3.2 normalisation fix), and l = the smallest length-scale to
# resolve. Hake true (lx, ly, lz) = (50, 300, 150) gives min M of
# (24, 10, 56). v3.2 keeps M_x / M_y at v3 values for now and bumps
# M_z from 8 to 20 - that's a ~88 m basis resolution against a 150 m
# z length-scale, comfortable without exploding the parameter count.
# (M_z = 56 is the textbook target; we step there if recovery still
# fights us on the bathymetric structure.)
HSGP_M <- c(10L, 8L, 20L)    # M_total = 1600  (was 640 in v3)
HSGP_C <- c(1.5, 1.5, 1.5)   # boundary extension (>= 1.5 recommended)

OUTPUT_DIR <- "outputs/whale_edna_output_v3.2"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load simulation
# -----------------------------------------------------------------------------
cat("=== Step 1: Loading simulation ===\n")
sim <- readRDS("outputs/whale_edna_output_v3.2/whale_edna_sim_v3.2.rds")

samples         <- sim$design$samples
stations        <- sim$design$stations
gp_params       <- sim$truth$gp_params
lambda_true_si  <- sim$truth$lambda_true_si
zsample_effect  <- sim$truth$zsample_effect
C_obs_si        <- sim$truth$C_obs_si
qpcr_detect_raw <- sim$observed$qpcr_detect_raw
qpcr_ct_raw     <- sim$observed$qpcr_ct_raw

N            <- sim$meta$N
S            <- sim$meta$n_species
R            <- sim$meta$n_qpcr_rep
sp_common    <- sim$meta$sp_common
vol_aliquot  <- 2        # microlitres (from simulation script)
conv_factor  <- 10       # animals -> copies

cat(sprintf("  Loaded: N=%d samples, S=%d species, R=%d qPCR reps\n", N, S, R))
stopifnot(S == 1L)

# -----------------------------------------------------------------------------
# 2. Water column log-effect (zero everywhere — surface only)
# -----------------------------------------------------------------------------
cat("=== Step 2: Preparing log_zsample_effect ===\n")

log_zsample_effect <- log(as.matrix(zsample_effect))   # N x S, all zeros

log_zsample_effect[is.nan(log_zsample_effect)]                               <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect < 0] <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect > 0] <-  10.0

cat(sprintf("  log_zsample_effect range: [%.2f, %.2f]\n",
            min(log_zsample_effect), max(log_zsample_effect)))
stopifnot(all(is.finite(log_zsample_effect)))

# -----------------------------------------------------------------------------
# 3. GP coordinates - normalise to [-1, 1]
# -----------------------------------------------------------------------------
cat("=== Step 3: Normalising coordinates ===\n")

coords_raw <- cbind(
  as.numeric(samples[["X"]]),
  as.numeric(samples[["Y"]]),
  as.numeric(samples[["Z_bathy"]])
)

# v3.2 fix: normalisation constants are derived from the actual v3
# domain (500 km x 1270 km x ~3500 m), not the v1/v2 300 x 400 km
# legacy values that were never updated when v3 extended the domain.
# With L_hsgp = 1.5 the HSGP basis is valid only for normalised
# coordinates in [-1.5, 1.5] - under the old constants Y_norm reached
# 5.35, putting most northern stations outside the basis support.
X_km_max    <- sim$meta$X_km_max                  # 500
Y_km_max    <- sim$meta$Y_km_max                  # 1270
Z_bathy_max <- 3500                                # >= max(Z_bathy) with headroom

coord_centre <- c(X_km_max  / 2, Y_km_max  / 2, Z_bathy_max / 2)
coord_scale  <- c(X_km_max  / 2, Y_km_max  / 2, Z_bathy_max / 2)

coords_norm <- sweep(coords_raw, 2, coord_centre, "-")
coords_norm <- sweep(coords_norm, 2, coord_scale,  "/")

cat(sprintf("  coord_centre: [%.0f, %.0f, %.0f]\n",
            coord_centre[1], coord_centre[2], coord_centre[3]))
cat(sprintf("  coord_scale:  [%.0f, %.0f, %.0f]\n",
            coord_scale[1],  coord_scale[2],  coord_scale[3]))
cat(sprintf("  X_norm:       [%.3f, %.3f]\n", min(coords_norm[,1]), max(coords_norm[,1])))
cat(sprintf("  Y_norm:       [%.3f, %.3f]\n", min(coords_norm[,2]), max(coords_norm[,2])))
cat(sprintf("  Z_bathy_norm: [%.3f, %.3f]\n", min(coords_norm[,3]), max(coords_norm[,3])))
stopifnot(all(is.finite(coords_norm)))
stopifnot(max(abs(coords_norm)) < 1.5)             # all data inside HSGP boundary

# -----------------------------------------------------------------------------
# 4. qPCR data - long form (N * R rows)
# -----------------------------------------------------------------------------
cat("=== Step 4: Preparing qPCR data ===\n")

N_qpcr_long     <- N * R
qpcr_sample_idx <- rep(seq_len(N), each = R)

qpcr_detect_vec <- as.integer(qpcr_detect_raw[, 1])
qpcr_ct_vec     <- as.numeric(qpcr_ct_raw[, 1])
qpcr_ct_vec[is.na(qpcr_ct_vec)] <- 0.0   # ignored when detect == 0

stopifnot(all(qpcr_detect_vec %in% c(0L, 1L)))
stopifnot(length(qpcr_detect_vec) == N_qpcr_long)

cat(sprintf("  N_qpcr_long: %d\n", N_qpcr_long))
cat(sprintf("  qPCR detection rate (hake): %.1f%%\n",
            100 * mean(qpcr_detect_vec)))

# -----------------------------------------------------------------------------
# 5. Assemble Stan data list
# -----------------------------------------------------------------------------
cat("=== Step 5: Assembling Stan data list ===\n")

M_total <- prod(HSGP_M)

stan_data <- list(

  # Dimensions
  N = N,
  S = S,
  M = M_total,

  # Coordinates (normalised to [-1, 1])
  coords      = coords_norm,
  coord_scale = coord_scale,

  # HSGP
  L_hsgp = HSGP_C,
  m_hsgp = HSGP_M,

  # Water column offset (fixed; all zeros at surface)
  log_zsample_effect = log_zsample_effect,

  # Conversion factor
  log_conv_factor = log(conv_factor),

  # Aliquot volume (microlitres)
  vol_aliquot = vol_aliquot,

  # qPCR (hake only, long form)
  N_qpcr_long     = N_qpcr_long,
  qpcr_sample_idx = qpcr_sample_idx,
  qpcr_detect     = qpcr_detect_vec,
  qpcr_ct         = qpcr_ct_vec,

  # Fixed qPCR calibration. v3.2: kappa joins alpha_ct/beta_ct as data.
  alpha_ct = sim$truth$qpcr_params$alpha_ct,
  beta_ct  = sim$truth$qpcr_params$beta_ct,
  kappa    = sim$truth$qpcr_params$kappa,

  # Prior hyperparameters (Normal mu / sigma for each scalar or vector prior)
  prior_mu_sp_mu         =   2.0,
  prior_mu_sp_sig        =   1.5,
  prior_gp_sigma_sig     =   1.5,
  prior_gp_lx_mu         =  50.0, prior_gp_lx_sig        =  40.0,
  prior_gp_ly_mu         = 150.0, prior_gp_ly_sig        =  80.0,
  prior_gp_lz_mu         = 300.0, prior_gp_lz_sig        = 150.0,
  prior_sigma_ct_mu      =   0.5, prior_sigma_ct_sig     =   0.3
)

# Sanity checks
stopifnot(nrow(coords_norm)        == N)
stopifnot(nrow(log_zsample_effect) == N)
stopifnot(ncol(log_zsample_effect) == S)
stopifnot(length(qpcr_sample_idx)  == N_qpcr_long)

cat(sprintf("  Stan data ready: N=%d  S=%d  M=%d  N_qpcr=%d\n",
            N, S, M_total, N_qpcr_long))

# -----------------------------------------------------------------------------
# 6. Save
# -----------------------------------------------------------------------------
saveRDS(stan_data, file.path(OUTPUT_DIR, "stan_data.rds"))
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "stan_data.rds")))
