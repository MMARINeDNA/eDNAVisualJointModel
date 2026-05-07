# =============================================================================
# 03_format_stan_data_v1.R
#
# Read the simulated eDNA data (outputs/whale_edna_sim_v1.rds) and assemble
# the Stan data list used by stan/whale_edna_hsgp_v1.stan. Writes:
#
#   outputs/whale_edna_output_v1/stan_data.rds
# =============================================================================

library(tidyverse)

set.seed(42)

# HSGP basis counts
HSGP_M <- c(6L, 7L, 6L)
HSGP_C <- c(1.3, 1.3, 1.3)

OUTPUT_DIR <- "outputs/whale_edna_output_v1"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Step 1: Loading simulation ===\n")
sim <- readRDS("outputs/whale_edna_sim_v1.rds")

samples       <- sim$design$samples
stations      <- sim$design$stations
N             <- sim$meta$N
S             <- sim$meta$n_species
R             <- sim$meta$n_qpcr_rep
K             <- sim$meta$n_mb_rep
sp_names      <- sim$meta$sp_names
sp_common     <- sim$meta$sp_common
sample_depths <- sim$meta$sample_depths
vol           <- rep(sim$meta$vol_filtered, N)

cat(sprintf("  Loaded: N=%d samples, S=%d species, R=%d qPCR reps, K=%d MB reps\n",
            N, S, R, K))

cat("=== Step 2: Preparing Stan data ===\n")

# GP coordinates: (X, Y, Z_bathy)
coords_raw <- as.matrix(samples[, c("X", "Y", "Z_bathy")])

# Normalise each dimension to [-1, 1] using fixed physical half-ranges
coord_center <- c(150,  200,  1250)
coord_scale  <- c(150,  200,  1250)
coords_norm  <- sweep(coords_raw, 2, coord_center, "-")
coords_norm  <- sweep(coords_norm, 2, coord_scale,  "/")

L_hsgp  <- HSGP_C
M_total <- prod(HSGP_M)

cat(sprintf("  Coordinate ranges after normalisation:\n"))
cat(sprintf("    X:       [%.2f, %.2f]\n", min(coords_norm[,1]), max(coords_norm[,1])))
cat(sprintf("    Y:       [%.2f, %.2f]\n", min(coords_norm[,2]), max(coords_norm[,2])))
cat(sprintf("    Z_bathy: [%.2f, %.2f]\n", min(coords_norm[,3]), max(coords_norm[,3])))

# Water column effect: log_zsample_effect
log_zsample_effect <- log(sim$truth$zsample_effect)

# qPCR data (hake only)
n_detect    <- sim$observed$qpcr_n_detect
any_detect  <- sim$observed$qpcr_any_detect
mean_ct_obs <- sim$observed$qpcr_mean_ct
mean_ct_obs[is.na(mean_ct_obs)] <- 0.0

# Metabarcoding
mb_reads_arr <- sim$observed$mb_reads_rep    # N × S × K
mb_total_arr <- sim$observed$mb_total_rep    # N × K

# Guards
vol <- pmax(vol, 1e-6)
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect < 0] <- -10.0
log_zsample_effect[is.infinite(log_zsample_effect) & log_zsample_effect > 0] <-  10.0
log_zsample_effect[is.nan(log_zsample_effect)] <- -10.0
cat("  vol range after guard:",               range(vol), "\n")
cat("  log_zsample_effect range after guard:", range(log_zsample_effect), "\n")

stan_data <- list(
  N                  = N,
  S                  = S,
  R                  = R,
  K                  = K,
  M                  = M_total,
  coords             = coords_norm,
  coord_scale        = coord_scale,
  L_hsgp             = L_hsgp,
  m_hsgp             = HSGP_M,
  vol                = vol,
  log_zsample_effect = log_zsample_effect,
  n_detect           = n_detect,
  mean_ct_obs        = mean_ct_obs,
  any_detect         = any_detect,
  mb_reads           = mb_reads_arr,
  mb_total           = mb_total_arr
)

cat(sprintf("  Stan data prepared: N=%d, S=%d, R=%d, K=%d, M=%d\n",
            N, S, R, K, M_total))

saveRDS(stan_data, file.path(OUTPUT_DIR, "stan_data.rds"))
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "stan_data.rds")))
