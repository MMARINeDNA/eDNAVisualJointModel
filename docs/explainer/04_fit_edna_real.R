# eDNA model on the real MURI bottles: the EARLIER fit, made before the qPCR lab's
# file arrived, when qPCR had no Ct or standard curve. It fits metabarcoding only:
# all qPCR rows are dropped except one non-detect, because hsgp_2d_bathysp.stan
# requires N_qpcr_long >= 1. The eDNA object now carries Ct, the standards,
# dilutions, volumes and target depths; this script is to be rewritten to use them
# (estimated per-plate curves, dilution offsets, per-bottle volume, depth smooth;
# ROADMAP_2026-10-21, section 5 item 7).
# Run from the repo root: Rscript docs/explainer/04_fit_edna_real.R
source("docs/explainer/00_common.R")
source(file.path(HSGP4EDNA_DIR, "R/functions.R"))
ed <- readRDS("data/processed/muri_edna.rds")
o <- ed$observed
# one non-detect on a bottle with no hake reads in any MB replicate
hake_by_bottle <- tapply(o$mb_reads[, "hake"], o$mb_sample_idx, sum)
zero_b <- as.integer(names(hake_by_bottle)[hake_by_bottle == 0])
keep <- which(o$qpcr_sample_idx %in% zero_b & o$qpcr_detect == 0)[1]
stopifnot(!is.na(keep))
o$qpcr_sample_idx <- o$qpcr_sample_idx[keep]; o$qpcr_detect <- o$qpcr_detect[keep]
o$qpcr_ct <- 0
ed$observed <- o
ed$meta$N_qpcr_long <- 1L
# qPCR calibration is unknown (NA in the object). The single dummy non-detect makes the
# values irrelevant in practice; use the simulators' values so the data list is complete.
ed$truth$qpcr_params <- list(kappa = 0.85, alpha_ct = 38.0, beta_ct = 1.44,
                             gamma_0 = 0.30, gamma_1 = -0.30, sigma_0 = 0.23)
fd <- format_stan_data_gp2d_bathysp(ed, gp_prior = "design", M_max = GP_BASIS_BUDGET,   # 800 from v1.4 (was 400)
                                    bottle_noise = "gamma")
sd <- fd$stan_data
S <- ed$meta$n_species; D1 <- sd$D1; Mtot <- sd$M; K <- sd$K_bathy
cat(sprintf("N=%d S=%d M=%d (%dx%d) K=%d N_qpcr=%d N_mb=%d\n", sd$N, S, Mtot, sd$m_hsgp[1],
            sd$m_hsgp[2], K, sd$N_qpcr_long, sd$N_mb_long))
saveRDS(sd, file.path(EXPL_OUT, "edna_real_standata.rds"))
mod <- edna_model(file.path(HSGP4EDNA_DIR, "stan/hsgp_2d_bathysp.stan"),   # recompiles when includes change
                  file.path(HSGP4EDNA_DIR, "stan/include"))
init_fn <- function() c(list(mu_sp = c(runif(S, 1.8, 2.0), runif(1, 4.0, 4.1)),
  gp_sigma = runif(S, 0.8, 1.2), gp_l_raw = matrix(runif(S * D1, 0.2, 0.4), S, D1, byrow = TRUE),
  z_beta = matrix(0, S, Mtot), beta_bathy = matrix(0, S, K),
  beta0_phi = 0, gamma0_phi = 2.0, gamma1_phi = 1.0), bottle_noise_init(sd))
set.seed(42L); t0 <- Sys.time()
fit <- mod$sample(data = sd, chains = 2, parallel_chains = 2, iter_warmup = 400, iter_sampling = 400,
                  adapt_delta = 0.85, max_treedepth = 12, seed = 42L, init = init_fn, refresh = 50)
cat("runtime (min):", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "\n")
fit$save_object(file.path(EXPL_OUT, "edna_real_fit.rds"))
print(fit$diagnostic_summary())
