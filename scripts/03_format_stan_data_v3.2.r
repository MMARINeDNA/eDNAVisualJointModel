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

# HSGP basis terms per dimension. Riutort-Mayol et al. (2023):
#   m >= 1.75 * c / rho_l  where rho_l = l / S_norm and S_norm = 1.
# With c = 1.5 and the v3.2 normalisation fix (half-ranges 250 km,
# 635 km, 1750 m), the named (lx, ly, lz) = (50, 300, 150) require
# m >= (13, 6, 31). v3.2 now sets M = (14, 8, 32) so every dimension
# sits at or above the faithful-representation threshold. This is paired
# with the Option-A simulation change (zero-mean GP, no preference
# function) so length-scale recovery is a clean test - the truth is
# unambiguously a single SE GP with these length-scales, and the basis
# is rich enough to represent it.
HSGP_M <- c(14L, 8L, 32L)    # M_total = 3584  (was 1600)
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
vol_filtered <- sim$meta$vol_filtered   # litres of seawater filtered (= 2.5 L)

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
  coords       = coords_norm,
  coord_centre = coord_centre,   # not used by Stan; needed for post-hoc grid evals
  coord_scale  = coord_scale,

  # HSGP
  L_hsgp = HSGP_C,
  m_hsgp = HSGP_M,

  # Water column offset (fixed; all zeros at surface)
  log_zsample_effect = log_zsample_effect,

  # Conversion factor + volume filtered (both contribute to expected
  # bottle copies = conv_factor * lambda * zsample_effect * vol_filtered).
  log_conv_factor   = log(conv_factor),
  log_vol_filtered  = log(vol_filtered),

  # Aliquot volume (microlitres)
  vol_aliquot = vol_aliquot,

  # qPCR (hake only, long form)
  N_qpcr_long     = N_qpcr_long,
  qpcr_sample_idx = qpcr_sample_idx,
  qpcr_detect     = qpcr_detect_vec,
  qpcr_ct         = qpcr_ct_vec,

  # Fixed qPCR calibration. v3.2: kappa and sigma_ct now join
  # alpha_ct/beta_ct as data - the entire standard-curve fit is treated
  # as known. sigma_ct deliberately set HIGHER than the sim's PCR-noise
  # truth (0.4): the simulated Ct's also carry discrete-count noise from
  # Ct ~ log(integer count), which contributes ~0.46 of additional SD
  # at the simulated mean aliquot copies (E[C] ~ 10, Var[log(C)] ~ 0.1,
  # times beta_ct = 1.44). Pinning sigma_ct = 0.4 (PCR-only) made the
  # model fit the field too tightly to absorb that extra noise -
  # 50% max-treedepth saturation, length-scales biased low. 0.6 is the
  # honest combined PCR + discrete-count noise (sqrt(0.4^2 + 0.46^2)),
  # giving the field a noise budget for what it can't explain.
  alpha_ct = sim$truth$qpcr_params$alpha_ct,
  beta_ct  = sim$truth$qpcr_params$beta_ct,
  kappa    = sim$truth$qpcr_params$kappa,
  sigma_ct = 0.6,

  # Prior hyperparameters. gp_sigma uses Gamma(shape, rate). Tightened
  # from Gamma(4, 2) (mode 1.5, mean 2.0, sd 1.0) to Gamma(8, 4) (mode
  # 1.75, mean 2.0, sd 0.71) - same mean, much less mass below ~0.3.
  # The looser prior was visibly bimodal across chains: some chains
  # found the correct gp_sigma ~ 1 mode, others got trapped near zero
  # (Rhat 1.74; median 0.008, q95 1.99). With M_total = 3584 the
  # posterior has a near-degenerate sigma <-> z_beta ridge connecting
  # the two modes through a funnel-shaped region NUTS can't traverse
  # during warmup. The tighter prior makes the low-sigma mode
  # essentially unreachable.
  prior_mu_sp_mu         =   2.0,
  prior_mu_sp_sig        =   1.5,
  prior_gp_sigma_shape   =   8.0,
  prior_gp_sigma_rate    =   4.0,
  # gp_l priors. The previous (50, 150, 300) means were leftovers from
  # an earlier sim configuration (likely v1/v2) and never updated when
  # v3 set (lx, ly, lz) = (50, 300, 150). With M_total = 3584 the
  # length-scale posterior is mildly multi-modal and even small prior
  # mis-centring was producing chains stuck at different gp_l values
  # (Rhat 1.4+ on lx and lz). Centred on the v3 truth and tightened
  # ~30% to fight the multi-modality.
  prior_gp_lx_mu         =  50.0, prior_gp_lx_sig        =  30.0,
  prior_gp_ly_mu         = 300.0, prior_gp_ly_sig        = 100.0,
  prior_gp_lz_mu         = 150.0, prior_gp_lz_sig        =  50.0
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
