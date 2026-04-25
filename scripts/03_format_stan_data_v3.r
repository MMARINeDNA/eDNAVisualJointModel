# =============================================================================
# 03_format_stan_data_v3.R
#
# Read the simulated eDNA data (outputs/whale_edna_sim_v3.rds) and assemble
# the Stan data list used by stan/whale_edna_hsgp_v3.stan. Writes:
#
#   outputs/whale_edna_output_v3/stan_data.rds
#
# Data preparation:
#   - GP coordinates: (X, Y, Z_bathy), normalised to [-1, 1]
#   - log_zsample_effect: fixed N × S water-column offset, pre-computed in R
#   - qPCR: long-form N*R rows (one per replicate), hake only
#   - MB: long-form N*K rows (one per replicate), filtered to rows where
#         at least one species has reads > 0
# =============================================================================

library(tidyverse)

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------

# HSGP basis terms per dimension.
# lx_min ~ 30 km in 300 km domain  → ~10 terms
# ly_min ~ 100 km in 400 km domain → ~8 terms
# lz_min ~ 100 m in 2500 m domain  → ~8 terms
HSGP_M <- c(10L, 8L, 8L)     # M_total = 640
HSGP_C <- c(1.5, 1.5, 1.5)   # boundary extension (>= 1.5 recommended)

OUTPUT_DIR <- "outputs/whale_edna_output_v3"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load simulation
# -----------------------------------------------------------------------------
cat("=== Step 1: Loading simulation ===\n")
sim <- readRDS("outputs/whale_edna_sim_v3.rds")

samples         <- sim$design$samples
stations        <- sim$design$stations
gp_params       <- sim$truth$gp_params
lambda_true_si  <- sim$truth$lambda_true_si
zsample_effect  <- sim$truth$zsample_effect
C_obs_si        <- sim$truth$C_obs_si
qpcr_detect_raw <- sim$observed$qpcr_detect_raw
qpcr_ct_raw     <- sim$observed$qpcr_ct_raw
mb_reads_rep    <- sim$observed$mb_reads_rep

N            <- sim$meta$N
S            <- sim$meta$n_species
R            <- sim$meta$n_qpcr_rep
K            <- sim$meta$n_mb_rep
sp_common    <- sim$meta$sp_common
vol_aliquot  <- 2        # microlitres (from simulation script)
conv_factor  <- 10       # animals → copies

cat(sprintf("  Loaded: N=%d samples, S=%d species, R=%d qPCR reps, K=%d MB reps\n",
            N, S, R, K))

# -----------------------------------------------------------------------------
# 2. Water column log-effect (fixed offset, pre-computed from simulation)
# -----------------------------------------------------------------------------
cat("=== Step 2: Preparing log_zsample_effect ===\n")

log_zsample_effect <- log(as.matrix(zsample_effect))   # N × S

# Guard against -Inf / NaN at extreme depths
log_zsample_effect[is.nan(log_zsample_effect)]                               <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect < 0] <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect > 0] <-  10.0

cat(sprintf("  log_zsample_effect range: [%.2f, %.2f]\n",
            min(log_zsample_effect), max(log_zsample_effect)))
stopifnot(all(is.finite(log_zsample_effect)))

# -----------------------------------------------------------------------------
# 3. GP coordinates — normalise to [-1, 1]
# -----------------------------------------------------------------------------
cat("=== Step 3: Normalising coordinates ===\n")

coords_raw <- cbind(
  as.numeric(samples[["X"]]),
  as.numeric(samples[["Y"]]),
  as.numeric(samples[["Z_bathy"]])
)

# Fixed physical ranges (interpretable; consistent across subsets)
# X:       0 – 300 km  → centre = 150, half-range = 150
# Y:       0 – 400 km  → centre = 200, half-range = 200
# Z_bathy: 0 – 2500 m  → centre = 1250, half-range = 1250
coord_centre <- c(150,  200,  1250)
coord_scale  <- c(150,  200,  1250)

coords_norm <- sweep(coords_raw, 2, coord_centre, "-")
coords_norm <- sweep(coords_norm, 2, coord_scale,  "/")

cat(sprintf("  X_norm:       [%.3f, %.3f]\n", min(coords_norm[,1]), max(coords_norm[,1])))
cat(sprintf("  Y_norm:       [%.3f, %.3f]\n", min(coords_norm[,2]), max(coords_norm[,2])))
cat(sprintf("  Z_bathy_norm: [%.3f, %.3f]\n", min(coords_norm[,3]), max(coords_norm[,3])))
stopifnot(all(is.finite(coords_norm)))

# -----------------------------------------------------------------------------
# 4. qPCR data — long form (N * R rows)
# -----------------------------------------------------------------------------
cat("=== Step 4: Preparing qPCR data ===\n")

N_qpcr_long     <- N * R
qpcr_sample_idx <- rep(seq_len(N), each = R)

qpcr_detect_vec <- as.integer(qpcr_detect_raw[, 1])
qpcr_ct_vec     <- as.numeric(qpcr_ct_raw[, 1])
qpcr_ct_vec[is.na(qpcr_ct_vec)] <- 0.0   # Stan requires numeric; ignored when detect == 0

stopifnot(all(qpcr_detect_vec %in% c(0L, 1L)))
stopifnot(length(qpcr_detect_vec) == N_qpcr_long)

cat(sprintf("  N_qpcr_long: %d\n", N_qpcr_long))
cat(sprintf("  qPCR detection rate (hake): %.1f%%\n",
            100 * mean(qpcr_detect_vec)))

# -----------------------------------------------------------------------------
# 5. Metabarcoding data — long form (N * K rows), zero-read rows filtered
# -----------------------------------------------------------------------------
cat("=== Step 5: Preparing metabarcoding data ===\n")

# mb_reads_rep from rmultinom is S × (N*K); transpose to (N*K) × S
mb_reads_long_full <- t(mb_reads_rep)
storage.mode(mb_reads_long_full) <- "integer"

mb_sample_idx_full <- rep(seq_len(N), each = K)
mb_total_full      <- rowSums(mb_reads_long_full)

# Filter rows where ALL species have zero reads
keep_mb      <- mb_total_full > 0
n_dropped_mb <- sum(!keep_mb)

mb_reads_long <- mb_reads_long_full[keep_mb, , drop = FALSE]
mb_sample_idx <- mb_sample_idx_full[keep_mb]
mb_total_vec  <- mb_total_full[keep_mb]
N_mb_long     <- nrow(mb_reads_long)
storage.mode(mb_total_vec) <- "integer"

cat(sprintf("  MB rows (N*K = %d): %d kept, %d dropped (zero reads)\n",
            N * K, N_mb_long, n_dropped_mb))
cat(sprintf("  MB read depth range: %d – %d\n",
            min(mb_total_vec), max(mb_total_vec)))

stopifnot(all(mb_total_vec > 0))
stopifnot(all(mb_reads_long >= 0L))
stopifnot(ncol(mb_reads_long) == S)

# -----------------------------------------------------------------------------
# 6. Assemble Stan data list
# -----------------------------------------------------------------------------
cat("=== Step 6: Assembling Stan data list ===\n")

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

  # Water column offset (fixed)
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

  # Metabarcoding (all species, long form, zero-read rows removed)
  N_mb_long     = N_mb_long,
  mb_sample_idx = mb_sample_idx,
  mb_reads      = mb_reads_long,
  mb_total      = mb_total_vec,

  # Zero-inflation switch
  use_zi = 1L,

  # Fixed qPCR standard-curve coefficients. Sourced from the simulation
  # truth here (so the fit uses the true calibration); in a real analysis
  # these would come from a pre-estimated standard curve.
  alpha_ct = sim$truth$qpcr_params$alpha_ct,
  beta_ct  = sim$truth$qpcr_params$beta_ct,

  # Prior hyperparameters (Normal mu / sigma for each scalar or vector prior)
  prior_mu_sp_mu         =   2.0,
  prior_mu_sp_sig        =   1.5,
  prior_gp_sigma_sig     =   1.5,
  prior_gp_lx_mu         =  50.0, prior_gp_lx_sig        =  40.0,
  prior_gp_ly_mu         = 150.0, prior_gp_ly_sig        =  80.0,
  prior_gp_lz_mu         = 300.0, prior_gp_lz_sig        = 150.0,
  prior_kappa_mu         =   0.5, prior_kappa_sig        =   0.3,
  prior_sigma_ct_mu      =   0.5, prior_sigma_ct_sig     =   0.3,
  prior_beta0_phi_mu     =   0.0, prior_beta0_phi_sig    =   1.0,
  prior_gamma0_phi_mu    =   2.0, prior_gamma0_phi_sig   =   1.0,
  prior_gamma1_phi_mu    =   0.5, prior_gamma1_phi_sig   =   0.3
)

# Sanity checks
stopifnot(nrow(coords_norm)        == N)
stopifnot(nrow(log_zsample_effect) == N)
stopifnot(ncol(log_zsample_effect) == S)
stopifnot(length(qpcr_sample_idx)  == N_qpcr_long)
stopifnot(length(mb_sample_idx)    == N_mb_long)
stopifnot(nrow(mb_reads_long)      == N_mb_long)

cat(sprintf("  Stan data ready: N=%d  S=%d  M=%d  N_qpcr=%d  N_mb=%d\n",
            N, S, M_total, N_qpcr_long, N_mb_long))

# -----------------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------------
saveRDS(stan_data, file.path(OUTPUT_DIR, "stan_data.rds"))
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "stan_data.rds")))
