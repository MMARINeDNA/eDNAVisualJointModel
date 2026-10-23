# Simulate bathysp_depthspref (40 stations) and fit the 2-D HSGP + spline eDNA model.
# Run from the repo root: Rscript docs/explainer/01_fit_edna_sim.R
source("docs/explainer/00_common.R")
source(file.path(HSGP4EDNA_DIR, "R/functions.R"))

sample_depths <- c(0, 50, 150, 250, 350, 500)
zsample_pref  <- list(hake     = c(0.0,  1.0,  3.0,  2.57, 2.14,  1.50),
                      humpback = c(0.0,  0.67, 2.00, 1.06, 0.11, -1.30),
                      pwsd     = c(0.0, -1.0, -3.0, -5.0, -7.0, -10.0))
sim <- simulate_bathysp(seed = 303L, n_stations = 40L, sample_depths = sample_depths,
                        zsample_pref = zsample_pref)
saveRDS(sim, file.path(EXPL_OUT, "edna_sim.rds"))

fd <- format_stan_data_gp2d_bathysp(sim, gp_prior = "design", M_max = 150L, bottle_noise = "gamma")
sd <- fd$stan_data
S <- sim$meta$n_species; D1 <- sd$D1; Mtot <- sd$M; K <- sd$K_bathy
cat(sprintf("N=%d S=%d M=%d (%dx%d) K=%d N_qpcr=%d N_mb=%d\n", sd$N, S, Mtot, sd$m_hsgp[1],
            sd$m_hsgp[2], K, sd$N_qpcr_long, sd$N_mb_long))
saveRDS(sd, file.path(EXPL_OUT, "edna_sim_standata.rds"))

mod <- edna_model(file.path(HSGP4EDNA_DIR, "stan/hsgp_2d_bathysp.stan"),   # recompiles when includes change
                  file.path(HSGP4EDNA_DIR, "stan/include"))
init_fn <- function() c(list(mu_sp = c(runif(S, 1.8, 2.0), runif(1, 4.0, 4.1)),
  gp_sigma = runif(S, 0.8, 1.2), gp_l_raw = matrix(runif(S * D1, 0.2, 0.4), S, D1, byrow = TRUE),
  z_beta = matrix(0, S, Mtot), beta_bathy = matrix(0, S, K),
  beta0_phi = 0, gamma0_phi = 2.0, gamma1_phi = 1.0), bottle_noise_init(sd))
set.seed(42L); t0 <- Sys.time()
fit <- mod$sample(data = sd, chains = 3, parallel_chains = 3, iter_warmup = 400, iter_sampling = 400,
                  adapt_delta = 0.85, max_treedepth = 10, seed = 42L, init = init_fn, refresh = 50)
cat("runtime (min):", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "\n")
fit$save_object(file.path(EXPL_OUT, "edna_sim_fit.rds"))
print(fit$diagnostic_summary())
