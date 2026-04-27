# =============================================================================
# 03_format_stan_data_v2.R
#
# Read the simulated eDNA data (outputs/whale_edna_output_v2/whale_edna_sim_v2.rds) and assemble
# the Stan data list used by stan/whale_edna_hsgp_v2.stan. Writes:
#
#   outputs/whale_edna_output_v2/stan_data.rds
# =============================================================================

library(tidyverse)

set.seed(42)

HSGP_M <- c(10L, 8L, 8L)
HSGP_C <- c(1.5, 1.5, 1.5)

OUTPUT_DIR <- "outputs/whale_edna_output_v2"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Step 1: Loading simulation ===\n")
sim <- readRDS("outputs/whale_edna_output_v2/whale_edna_sim_v2.rds")

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
vol_aliquot  <- 2
conv_factor  <- 10

cat(sprintf("  Loaded: N=%d samples, S=%d species, R=%d qPCR reps, K=%d MB reps\n",
            N, S, R, K))

cat("=== Step 2: Preparing log_zsample_effect ===\n")
log_zsample_effect <- log(as.matrix(zsample_effect))
log_zsample_effect[is.nan(log_zsample_effect)]                               <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect < 0] <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect > 0] <-  10.0
cat(sprintf("  log_zsample_effect range: [%.2f, %.2f]\n",
            min(log_zsample_effect), max(log_zsample_effect)))
stopifnot(all(is.finite(log_zsample_effect)))

cat("=== Step 3: Normalising coordinates ===\n")
coords_raw <- cbind(
  as.numeric(samples[["X"]]),
  as.numeric(samples[["Y"]]),
  as.numeric(samples[["Z_bathy"]])
)
coord_centre <- c(150,  200,  1250)
coord_scale  <- c(150,  200,  1250)
coords_norm <- sweep(coords_raw, 2, coord_centre, "-")
coords_norm <- sweep(coords_norm, 2, coord_scale,  "/")
cat(sprintf("  X_norm:       [%.3f, %.3f]\n", min(coords_norm[,1]), max(coords_norm[,1])))
cat(sprintf("  Y_norm:       [%.3f, %.3f]\n", min(coords_norm[,2]), max(coords_norm[,2])))
cat(sprintf("  Z_bathy_norm: [%.3f, %.3f]\n", min(coords_norm[,3]), max(coords_norm[,3])))
stopifnot(all(is.finite(coords_norm)))

cat("=== Step 4: Preparing qPCR data ===\n")
N_qpcr_long     <- N * R
qpcr_sample_idx <- rep(seq_len(N), each = R)
qpcr_detect_vec <- as.integer(qpcr_detect_raw[, 1])
qpcr_ct_vec     <- as.numeric(qpcr_ct_raw[, 1])
qpcr_ct_vec[is.na(qpcr_ct_vec)] <- 0.0
stopifnot(all(qpcr_detect_vec %in% c(0L, 1L)))
stopifnot(length(qpcr_detect_vec) == N_qpcr_long)
cat(sprintf("  N_qpcr_long: %d\n", N_qpcr_long))
cat(sprintf("  qPCR detection rate (hake): %.1f%%\n",
            100 * mean(qpcr_detect_vec)))

cat("=== Step 5: Preparing metabarcoding data ===\n")
mb_reads_long_full <- t(mb_reads_rep)
storage.mode(mb_reads_long_full) <- "integer"
mb_sample_idx_full <- rep(seq_len(N), each = K)
mb_total_full      <- rowSums(mb_reads_long_full)
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

cat("=== Step 6: Assembling Stan data list ===\n")
M_total <- prod(HSGP_M)

stan_data <- list(
  N = N, S = S, M = M_total,
  coords      = coords_norm,
  coord_scale = coord_scale,
  L_hsgp = HSGP_C,
  m_hsgp = HSGP_M,
  log_zsample_effect = log_zsample_effect,
  log_conv_factor = log(conv_factor),
  vol_aliquot = vol_aliquot,
  N_qpcr_long     = N_qpcr_long,
  qpcr_sample_idx = qpcr_sample_idx,
  qpcr_detect     = qpcr_detect_vec,
  qpcr_ct         = qpcr_ct_vec,
  N_mb_long     = N_mb_long,
  mb_sample_idx = mb_sample_idx,
  mb_reads      = mb_reads_long,
  mb_total      = mb_total_vec,
  use_zi = 1L,
  prior_mu_sp_mu     =  2.0,
  prior_mu_sp_sig    =  1.5,
  prior_gp_sigma_sig =  1.5,
  prior_gp_lx_mu     = 50.0,  prior_gp_lx_sig  = 40.0,
  prior_gp_ly_mu     = 150.0, prior_gp_ly_sig  = 80.0,
  prior_gp_lz_mu     = 300.0, prior_gp_lz_sig  = 150.0
)

stopifnot(nrow(coords_norm)        == N)
stopifnot(nrow(log_zsample_effect) == N)
stopifnot(ncol(log_zsample_effect) == S)
stopifnot(length(qpcr_sample_idx)  == N_qpcr_long)
stopifnot(length(mb_sample_idx)    == N_mb_long)
stopifnot(nrow(mb_reads_long)      == N_mb_long)

cat(sprintf("  Stan data ready: N=%d  S=%d  M=%d  N_qpcr=%d  N_mb=%d\n",
            N, S, M_total, N_qpcr_long, N_mb_long))

saveRDS(stan_data, file.path(OUTPUT_DIR, "stan_data.rds"))
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "stan_data.rds")))
