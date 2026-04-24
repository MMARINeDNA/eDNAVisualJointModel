# =============================================================================
# 03_format_stan_data_v4.R
#
# Read the simulated eDNA data (outputs/whale_edna_sim_v4.rds) and assemble
# the Stan data list used by stan/whale_edna_hsgp_v4.stan. Writes:
#
#   outputs/whale_edna_output_v4/stan_data.rds
#
# Compared to v3, the v4 sim emits the qPCR and metabarcoding data already
# in long form (one row per aliquot), with per-sample replicate counts that
# can vary between 1 and 3 (~20% of samples have fewer than 3 reps for each
# process). This script just threads those long-form vectors through to the
# Stan data list — the Stan model's long-form structure already handles
# variable replication via qpcr_sample_idx / mb_sample_idx.
# =============================================================================

library(tidyverse)

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
HSGP_M <- c(10L, 8L, 8L)     # M_total = 640
HSGP_C <- c(1.5, 1.5, 1.5)   # boundary extension (>= 1.5 recommended)

OUTPUT_DIR <- "outputs/whale_edna_output_v4"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load simulation
# -----------------------------------------------------------------------------
cat("=== Step 1: Loading simulation ===\n")
sim <- readRDS("outputs/whale_edna_sim_v4.rds")

samples        <- sim$design$samples
stations       <- sim$design$stations
gp_params      <- sim$truth$gp_params
lambda_true_si <- sim$truth$lambda_true_si
zsample_effect <- sim$truth$zsample_effect

N            <- sim$meta$N
S            <- sim$meta$n_species
sp_common    <- sim$meta$sp_common
vol_aliquot  <- sim$meta$vol_aliquot
conv_factor  <- sim$meta$conv_factor   # length-S vector, species-specific
stopifnot(length(conv_factor) == S)

N_qpcr_long  <- sim$meta$N_qpcr_long
N_mb_long    <- sim$meta$N_mb_long

cat(sprintf("  Loaded: N=%d samples, S=%d species\n", N, S))
cat(sprintf("  qPCR rep distribution: %s\n",
            paste(sprintf("%d:%d", as.integer(names(table(sim$design$n_qpcr_rep_i))),
                          as.integer(table(sim$design$n_qpcr_rep_i))),
                  collapse = ", ")))
cat(sprintf("  MB   rep distribution: %s\n",
            paste(sprintf("%d:%d", as.integer(names(table(sim$design$n_mb_rep_i))),
                          as.integer(table(sim$design$n_mb_rep_i))),
                  collapse = ", ")))

# -----------------------------------------------------------------------------
# 2. Water column log-effect
# -----------------------------------------------------------------------------
cat("=== Step 2: Preparing log_zsample_effect ===\n")

log_zsample_effect <- log(as.matrix(zsample_effect))   # N × S
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
coord_centre <- c(150,  200,  1250)
coord_scale  <- c(150,  200,  1250)
coords_norm  <- sweep(coords_raw, 2, coord_centre, "-")
coords_norm  <- sweep(coords_norm, 2, coord_scale,  "/")

cat(sprintf("  X_norm:       [%.3f, %.3f]\n", min(coords_norm[,1]), max(coords_norm[,1])))
cat(sprintf("  Y_norm:       [%.3f, %.3f]\n", min(coords_norm[,2]), max(coords_norm[,2])))
cat(sprintf("  Z_bathy_norm: [%.3f, %.3f]\n", min(coords_norm[,3]), max(coords_norm[,3])))
stopifnot(all(is.finite(coords_norm)))

# -----------------------------------------------------------------------------
# 4. qPCR long form — pass through
# -----------------------------------------------------------------------------
cat("=== Step 4: qPCR long form ===\n")

qpcr_sample_idx <- as.integer(sim$observed$qpcr_sample_idx)
qpcr_detect_vec <- as.integer(sim$observed$qpcr_detect)
qpcr_ct_vec     <- as.numeric(sim$observed$qpcr_ct)
qpcr_ct_vec[is.na(qpcr_ct_vec)] <- 0.0   # Stan requires numeric; ignored when detect == 0

stopifnot(length(qpcr_detect_vec) == N_qpcr_long)
stopifnot(all(qpcr_detect_vec %in% c(0L, 1L)))
stopifnot(all(qpcr_sample_idx >= 1L & qpcr_sample_idx <= N))

cat(sprintf("  N_qpcr_long: %d\n", N_qpcr_long))
cat(sprintf("  qPCR detection rate (hake): %.1f%%\n",
            100 * mean(qpcr_detect_vec)))

# -----------------------------------------------------------------------------
# 5. Metabarcoding long form — pass through, drop zero-read rows
# -----------------------------------------------------------------------------
cat("=== Step 5: Metabarcoding long form ===\n")

mb_reads_full <- sim$observed$mb_reads
storage.mode(mb_reads_full) <- "integer"
mb_sample_idx_full <- as.integer(sim$observed$mb_sample_idx)
mb_total_full      <- as.integer(sim$observed$mb_total)

keep_mb      <- mb_total_full > 0
n_dropped_mb <- sum(!keep_mb)

mb_reads_long <- mb_reads_full[keep_mb, , drop = FALSE]
mb_sample_idx <- mb_sample_idx_full[keep_mb]
mb_total_vec  <- mb_total_full[keep_mb]
N_mb_long     <- nrow(mb_reads_long)   # may have shrunk after dropping zeros
storage.mode(mb_total_vec) <- "integer"

cat(sprintf("  MB rows: %d kept, %d dropped (zero reads)\n",
            N_mb_long, n_dropped_mb))
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

  # Species-specific conversion factor (log scale), length S
  log_conv_factor = as.numeric(log(conv_factor)),

  # Aliquot volume (microlitres)
  vol_aliquot = vol_aliquot,

  # qPCR (hake only, long form — variable reps per sample)
  N_qpcr_long     = N_qpcr_long,
  qpcr_sample_idx = qpcr_sample_idx,
  qpcr_detect     = qpcr_detect_vec,
  qpcr_ct         = qpcr_ct_vec,

  # Metabarcoding (all species, long form — variable reps per sample)
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
