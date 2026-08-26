# =============================================================================
# demo_B_hsgp2D_bathycov_vS2.R
#
# Demonstration B: a 2-D HSGP over (X, Y) with bottom depth as a PARAMETRIC
# covariate (centred natural spline), instead of a 3rd GP axis. The truth
# (simulate_whale_edna_vS2B) is matched: log-density = mu + GP(X,Y) +
# B(Z_bathy) . beta. Demonstrates the covariate approach recovers both the
# spatial field and the smooth depth curve at a fraction of the 3-D cost.
#
# Outputs: outputs/whale_edna_output_vS2B/hsgp_bathycov_fit.rds
#          outputs/whale_edna_output_vS2B/hsgp_bathycov_recovery.pdf
#
# Env overrides: DEMO_MX, DEMO_MY, DEMO_CHAINS, DEMO_WARMUP, DEMO_SAMPLE, DEMO_SEED.
# =============================================================================

library(cmdstanr)
library(posterior)
suppressMessages(library(tidyverse))
source("scripts/simulation/vS1_functions.R")

OUTPUT_DIR <- "outputs/whale_edna_output_vS2B"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

MX <- as.integer(Sys.getenv("DEMO_MX", "14"))
MY <- as.integer(Sys.getenv("DEMO_MY", "6"))
N_CHAINS <- as.integer(Sys.getenv("DEMO_CHAINS", "3"))
N_WARMUP <- as.integer(Sys.getenv("DEMO_WARMUP", "300"))
N_SAMPLE <- as.integer(Sys.getenv("DEMO_SAMPLE", "400"))
SEED     <- as.integer(Sys.getenv("DEMO_SEED", "202"))

cat(sprintf("=== Demo B: 2-D HSGP (X,Y) + parametric spline Z_bathy ===\n"))
cat(sprintf("Basis: %d x %d = %d\n", MX, MY, MX*MY))

sim <- simulate_whale_edna_vS2B(seed = SEED)
saveRDS(sim, file.path(OUTPUT_DIR, "whale_edna_sim_vS2B.rds"))
S <- sim$meta$n_species
fd <- format_stan_data_bathycov(sim, HSGP_M = c(MX, MY))
sd <- fd$stan_data
D1 <- sd$D1; Mtot <- sd$M; K <- sd$K_bathy
cat(sprintf("N=%d  S=%d  M=%d  K_bathy=%d  N_mb=%d\n", sd$N, S, Mtot, K, sd$N_mb_long))

mod <- cmdstan_model("stan/whale_edna_hsgp_2D_bathycov.stan")
init_fn <- function() list(
  mu_sp = c(runif(S,1.8,2.0), runif(1,4.0,4.1)),
  gp_sigma = runif(S,0.8,1.2),
  gp_l_raw = matrix(runif(S*D1,0.2,0.4), S, D1, byrow=TRUE),
  z_beta = matrix(0, S, Mtot),
  beta_bathy = matrix(0, S, K),
  beta0_phi = 0, gamma0_phi = 2.0, gamma1_phi = 1.0)

fit <- mod$sample(data = sd, chains = N_CHAINS, parallel_chains = N_CHAINS,
                  iter_warmup = N_WARMUP, iter_sampling = N_SAMPLE,
                  adapt_delta = 0.80, max_treedepth = 11,
                  seed = SEED, init = init_fn, refresh = 50)

# -----------------------------------------------------------------------------
# Recovery: GP hyperparameters, spline coefficients, depth curve, field
# -----------------------------------------------------------------------------
gp <- sim$truth$gp_params
qsum <- summarise_draws(fit$draws(variables = c("gp_l","gp_sigma","mu_sp","beta_bathy")),
                        mean = mean,
                        q025 = ~as.numeric(stats::quantile(.x, 0.025)),
                        q975 = ~as.numeric(stats::quantile(.x, 0.975)), rhat = rhat)
getm <- function(v) qsum$mean[qsum$variable == v]

hyper <- map_dfr(seq_len(S), function(s) tibble(
  species = s, param = c("lx","ly","gp_sigma","mu"),
  true = c(gp[[s]]$lx, gp[[s]]$ly, gp[[s]]$sigma, gp[[s]]$mu),
  variable = c(sprintf("gp_l[%d,1]",s), sprintf("gp_l[%d,2]",s),
               sprintf("gp_sigma[%d]",s), sprintf("mu_sp[%d]",s)))) %>%
  left_join(dplyr::select(qsum, variable, mean, q025, q975), by="variable") %>%
  mutate(covered = true>=q025 & true<=q975) %>% dplyr::select(-variable)

# beta_bathy recovery
beta_true <- sim$truth$beta_bathy_true
beta_rec <- map_dfr(seq_len(S), function(s) map_dfr(seq_len(K), function(k) tibble(
  species=s, k=k, true=beta_true[s,k], mean=getm(sprintf("beta_bathy[%d,%d]",s,k)))))

# depth-curve recovery on a Z grid (using the sim's stored spline setup)
zgrid <- seq(min(sim$design$samples$Z_bathy), max(sim$design$samples$Z_bathy), length.out=100)
Bg <- bathy_spline_basis(zgrid, setup = sim$truth$bathy_spline)$B
beta_est <- matrix(0, S, K); for (s in 1:S) for (k in 1:K) beta_est[s,k] <- getm(sprintf("beta_bathy[%d,%d]",s,k))
curve_df <- map_dfr(seq_len(S), function(s) tibble(
  species = sim$meta$sp_common[s], Z_bathy = zgrid,
  true = as.numeric(Bg %*% beta_true[s,]), est = as.numeric(Bg %*% beta_est[s,])))

# latent field (GP + bathy)
ll <- summarise_draws(fit$draws(variables="log_lambda"), mean = mean)
m  <- stringr::str_match(ll$variable, "\\[(\\d+),(\\d+)\\]")
post <- matrix(NA_real_, sd$N, S); post[cbind(as.integer(m[,2]), as.integer(m[,3]))] <- ll$mean
truef <- log(sim$truth$lambda_true_si)
field <- map_dfr(seq_len(S), function(s) tibble(
  species=sim$meta$sp_common[s], R2=cor(post[,s],truef[,s])^2,
  rmse=sqrt(mean((post[,s]-truef[,s])^2)), bias=mean(post[,s]-truef[,s])))

diag <- fit$diagnostic_summary()
cat(sprintf("\nDivergences=%d  treedepth hits=%d  E-BFMI=%s\n", sum(diag$num_divergent),
            sum(diag$num_max_treedepth), paste(sprintf("%.2f",diag$ebfmi), collapse=", ")))
cat("\n--- GP hyperparameter recovery ---\n"); print(as.data.frame(hyper), digits=3, row.names=FALSE)
cat("\n--- beta_bathy recovery (spline coefficients) ---\n"); print(as.data.frame(beta_rec), digits=3, row.names=FALSE)
cat("\n--- Latent-field recovery ---\n"); print(as.data.frame(field), digits=3, row.names=FALSE)

saveRDS(list(hyper=hyper, beta_rec=beta_rec, curve_df=curve_df, field=field,
             diagnostic_summary=diag, HSGP_M=c(MX,MY)),
        file.path(OUTPUT_DIR, "hsgp_bathycov_fit.rds"))

pdf(file.path(OUTPUT_DIR, "hsgp_bathycov_recovery.pdf"), width=10, height=6, onefile=TRUE)
print(ggplot(curve_df, aes(Z_bathy)) +
        geom_line(aes(y=true, colour="true"), linewidth=1) +
        geom_line(aes(y=est, colour="estimated"), linewidth=1, linetype=2) +
        facet_wrap(~species) + scale_colour_manual(values=c(true="black", estimated="red")) +
        theme_bw() + labs(title="Demo B: bottom-depth effect recovery (parametric spline)",
                          x="Z_bathy (m)", y="h(Z_bathy) log-density effect", colour=NULL))
fdf <- map_dfr(seq_len(S), function(s) tibble(species=sim$meta$sp_common[s], true=truef[,s], post=post[,s]))
print(ggplot(fdf, aes(true, post)) + geom_point(alpha=0.5) +
        geom_abline(slope=1,intercept=0,colour="red") + facet_wrap(~species, scales="free") +
        theme_bw() + labs(title="Demo B: latent-field recovery (2-D HSGP + spline bathy)",
                          x="true log_lambda", y="posterior mean log_lambda"))
invisible(dev.off())
cat(sprintf("\nSaved %s and hsgp_bathycov_recovery.pdf\n", file.path(OUTPUT_DIR, "hsgp_bathycov_fit.rds")))
