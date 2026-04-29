# =============================================================================
# 03_format_stan_data_v4.2.R
#
# Read outputs/whale_edna_output_v4.2/whale_edna_sim_v4.2.rds and assemble
# the Stan data list used by stan/whale_edna_hsgp_v4.2.stan. Writes:
#
#   outputs/whale_edna_output_v4.2/stan_data.rds
#
# Changes from v4.1 

#   * HSGP_M sized to Ole's desires as opposed toRiutort-Mayol's faithful-representation rule for
#     the simulated truth: M = (10, 6, 16) covering lx=50 km, ly=300 km,
#     lz=150 m on the corrected coordinate scale.
#   * coord_centre / coord_scale derived from the sim's actual domain
#     extents (250 km, 635 km, 1750 m half-ranges) so all normalised
#     coords land in [-1, 1].
#   * gp_l priors centred on the sim truth: lx ~ N(50, 30),
#     ly ~ N(300, 100), lz ~ N(150, 50).
#   * gp_sigma prior is now Gamma(shape=8, rate=4): mode 1.75, mean 2.0,
#     sd 0.71. Replaces the half_normal that trapped gp_sigma at 0.
#   * log_vol_filtered passed as data and added to the eDNA log-mean
#     (was missing from v3 / v4).
#   * all qPCR parameters passed as data and perfectly match the 
#     simulated values.
# =============================================================================

library(tidyverse)

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
# HSGP basis terms per dimension. Sized to Riutort-Mayol's
# faithful-representation rule m >= 1.75 * c / rho_l for the simulated
# truth. With c = 1.5, half-ranges (250, 635, 1750), and named
# (lx, ly, lz) = (50, 300, 150) the rule wants m >= (13, 6, 31).
#HSGP_M <- c(14L, 8L, 32L)    # M_total = 3584
HSGP_M <- c(10L, 6L, 16L)    # M_total = 3584
HSGP_C <- c(1.5, 1.5, 1.5)   # boundary extension (>= 1.5 recommended)

OUTPUT_DIR <- "outputs/whale_edna_output_v4.2"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load simulation
# -----------------------------------------------------------------------------
cat("=== Step 1: Loading simulation ===\n")
sim <- readRDS("outputs/whale_edna_output_v4.2/whale_edna_sim_v4.2.rds")

samples        <- sim$design$samples
stations       <- sim$design$stations
gp_params      <- sim$truth$gp_params
lambda_true_si <- sim$truth$lambda_true_si
zsample_effect <- sim$truth$zsample_effect

N            <- sim$meta$N
S            <- sim$meta$n_species
sp_common    <- sim$meta$sp_common
vol_aliquot  <- sim$meta$vol_aliquot
vol_filtered <- sim$meta$vol_filtered                  # litres seawater filtered
conv_factor  <- sim$meta$conv_factor                   # length-S vector
stopifnot(length(conv_factor) == S+1)

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

# coord_centre / coord_scale derived from the actual sim domain extents
# so all normalised coords sit in [-1, 1] (and inside the HSGP boundary
# at L = 1.5).
X_km_max     <- sim$meta$X_km_max                  # 500
Y_km_max     <- sim$meta$Y_km_max                  # 1270
Z_bathy_max  <- 3500                                # >= max(Z_bathy) with headroom

coord_centre <- c(X_km_max  / 2, Y_km_max  / 2, Z_bathy_max / 2)
coord_scale  <- c(X_km_max  / 2, Y_km_max  / 2, Z_bathy_max / 2)

coords_norm  <- sweep(coords_raw, 2, coord_centre, "-")
coords_norm  <- sweep(coords_norm, 2, coord_scale,  "/")

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
stopifnot(ncol(mb_reads_long) == (S +1)) # +1 because 1 column for junk

# -----------------------------------------------------------------------------
# 6. Assemble Stan data list
# -----------------------------------------------------------------------------
# v4.2 now includes a junk observation column for the mb_reads.
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

  # log(volume filtered, litres). Now part of the eDNA log-mean
  # (was missing from v3 / v4).
  log_vol_filtered = log(vol_filtered),

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

  # Fixed qPCR standard-curve coefficients. v4.2: 
  # kappa, gamma_0, gamma_1, sigma_0 added as data.
  # this is the variance funnel for qPCR observations.
  alpha_ct = sim$truth$qpcr_params$alpha_ct,
  beta_ct  = sim$truth$qpcr_params$beta_ct,
  kappa    = sim$truth$qpcr_params$kappa,
  gamma0_ct    = sim$truth$qpcr_params$gamma_0,
  gamma1_ct    = sim$truth$qpcr_params$gamma_1,
  sigma0_ct    = sim$truth$qpcr_params$sigma_0,

  # Prior hyperparameters. gp_sigma now uses Gamma(shape, rate);
  # gp_l priors centred on the v4.2 truth (lx, ly, lz) = (50, 300, 150)
  # with ~30% SD. Other priors unchanged from v4.
  prior_mu_sp_mu         =   2.0,
  prior_mu_sp_sig        =   1.5,
  prior_gp_sigma_shape   =   8.0,
  prior_gp_sigma_rate    =   4.0,
  prior_gp_lx_mu         =  50.0, prior_gp_lx_sig        =  30.0,
  prior_gp_ly_mu         = 300.0, prior_gp_ly_sig        = 100.0,
  prior_gp_lz_mu         = 150.0, prior_gp_lz_sig        =  50.0,
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
