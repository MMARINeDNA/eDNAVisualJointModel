# =============================================================================
# demo_A_hsgp3D_vS2.R
#
# Demonstration A: the HSGP works in THREE dimensions, with bottom depth
# (Z_bathy) as the third GP axis. Fits the dimension-general model
# stan/whale_edna_hsgp_vS1.stan (D1 = 3) to the vS2 simulation (a genuine
# 3-D GP truth over X, Y, Z_bathy). The Z_bathy basis is CAPPED (default 20):
# the faithful-representation rule wants ~40 there (short bottom-depth length
# scale over a large range), which is impractical, so we accept some depth
# oversmoothing to keep total M tractable.
#
# Inputs : outputs/whale_edna_output_vS2/whale_edna_sim_vS2.rds  (run 01 vS2)
# Outputs: outputs/whale_edna_output_vS2/hsgp3D_fit.rds
#          outputs/whale_edna_output_vS2/hsgp3D_recovery.pdf
#
# Env overrides for smoke runs: DEMO_MX, DEMO_MY, DEMO_MZ, DEMO_CHAINS,
# DEMO_WARMUP, DEMO_SAMPLE.
# =============================================================================

library(cmdstanr)
library(posterior)
suppressMessages(library(tidyverse))
source("scripts/simulation/vS1_functions.R")

OUTPUT_DIR <- "outputs/whale_edna_output_vS2"
SIM_RDS    <- file.path(OUTPUT_DIR, "whale_edna_sim_vS2.rds")

# Basis: rule-ish in X, Y; Z_bathy capped for tractability.
MX <- as.integer(Sys.getenv("DEMO_MX", "14"))
MY <- as.integer(Sys.getenv("DEMO_MY", "6"))
MZ <- as.integer(Sys.getenv("DEMO_MZ", "20"))
N_CHAINS <- as.integer(Sys.getenv("DEMO_CHAINS", "3"))
N_WARMUP <- as.integer(Sys.getenv("DEMO_WARMUP", "300"))
N_SAMPLE <- as.integer(Sys.getenv("DEMO_SAMPLE", "400"))

cat(sprintf("=== Demo A: 3-D HSGP (X, Y, Z_bathy) on vS2 ===\n"))
cat(sprintf("Basis: %d x %d x %d = %d\n", MX, MY, MZ, MX*MY*MZ))

if (!file.exists(SIM_RDS))
  stop("Run scripts/simulation/01_simulate_whale_edna_vS2.r first (writes ", SIM_RDS, ")")
sim <- readRDS(SIM_RDS)
S   <- sim$meta$n_species

fd <- format_stan_data_hsgp(sim, HSGP_M = c(MX, MY, MZ),
                            coord_cols = c("X","Y","Z_bathy"))
sd <- fd$stan_data
D1 <- sd$D1; Mtot <- sd$M
cat(sprintf("N=%d  S=%d  M=%d  N_mb=%d  coord_scale=(%s)\n",
            sd$N, S, Mtot, sd$N_mb_long,
            paste(round(sd$coord_scale), collapse=", ")))

# -----------------------------------------------------------------------------
mod <- cmdstan_model("stan/whale_edna_hsgp_vS1.stan")
init_fn <- function() list(
  mu_sp = c(runif(S,1.8,2.0), runif(1,4.0,4.1)),
  gp_sigma = runif(S,0.8,1.2),
  gp_l_raw = matrix(runif(S*D1,0.2,0.4), S, D1, byrow=TRUE),
  z_beta = matrix(0, S, Mtot),
  beta0_phi = 0, gamma0_phi = 2.0, gamma1_phi = 1.0)

fit <- mod$sample(data = sd, chains = N_CHAINS, parallel_chains = N_CHAINS,
                  iter_warmup = N_WARMUP, iter_sampling = N_SAMPLE,
                  adapt_delta = 0.80, max_treedepth = 11,
                  seed = 42, init = init_fn, refresh = 50)

# -----------------------------------------------------------------------------
# Recovery: hyperparameters (lx, ly, lz, gp_sigma, mu) + latent field
# -----------------------------------------------------------------------------
gp <- sim$truth$gp_params
qsum <- summarise_draws(fit$draws(variables = c("gp_l","gp_sigma","mu_sp")),
                        mean = mean,
                        q025 = ~as.numeric(stats::quantile(.x, 0.025)),
                        q975 = ~as.numeric(stats::quantile(.x, 0.975)),
                        rhat = rhat)
truth_tbl <- map_dfr(seq_len(S), function(s) tibble(
  species = s, param = c("lx","ly","lz","gp_sigma","mu"),
  true = c(gp[[s]]$lx, gp[[s]]$ly, gp[[s]]$lz, gp[[s]]$sigma, gp[[s]]$mu),
  variable = c(sprintf("gp_l[%d,1]",s), sprintf("gp_l[%d,2]",s), sprintf("gp_l[%d,3]",s),
               sprintf("gp_sigma[%d]",s), sprintf("mu_sp[%d]",s))))
recovery <- truth_tbl %>%
  left_join(dplyr::select(qsum, variable, mean, q025, q975, rhat), by = "variable") %>%
  mutate(covered = true >= q025 & true <= q975) %>%
  dplyr::select(species, param, true, mean, q025, q975, covered, rhat)

# latent field
ll <- summarise_draws(fit$draws(variables = "log_lambda"), mean = mean)
m  <- stringr::str_match(ll$variable, "\\[(\\d+),(\\d+)\\]")
post <- matrix(NA_real_, sd$N, S); post[cbind(as.integer(m[,2]), as.integer(m[,3]))] <- ll$mean
truef <- log(sim$truth$lambda_true_si)
field <- map_dfr(seq_len(S), function(s) tibble(
  species = sim$meta$sp_common[s],
  R2 = cor(post[,s], truef[,s])^2,
  rmse = sqrt(mean((post[,s]-truef[,s])^2)), bias = mean(post[,s]-truef[,s])))

diag <- fit$diagnostic_summary()
cat(sprintf("\nDivergences=%d  max treedepth hits=%d  E-BFMI=%s\n",
            sum(diag$num_divergent), sum(diag$num_max_treedepth),
            paste(sprintf("%.2f", diag$ebfmi), collapse=", ")))
cat("\n--- Hyperparameter recovery (3-D: lx, ly, lz km/m) ---\n")
print(as.data.frame(recovery), digits = 3, row.names = FALSE)
cat("\n--- Latent-field recovery (log_lambda) ---\n")
print(as.data.frame(field), digits = 3, row.names = FALSE)

# -----------------------------------------------------------------------------
saveRDS(list(fit_draws = fit$draws(variables = c("lp__","mu_sp","gp_sigma","gp_l")),
             recovery = recovery, field = field, HSGP_M = c(MX,MY,MZ),
             diagnostic_summary = diag),
        file.path(OUTPUT_DIR, "hsgp3D_fit.rds"))

pdf(file.path(OUTPUT_DIR, "hsgp3D_recovery.pdf"), width = 10, height = 6, onefile = TRUE)
# field recovery scatter
fdf <- map_dfr(seq_len(S), function(s) tibble(
  species = sim$meta$sp_common[s], true = truef[,s], post = post[,s]))
print(ggplot(fdf, aes(true, post)) + geom_point(alpha=0.5) +
        geom_abline(slope=1, intercept=0, colour="red") +
        facet_wrap(~species, scales="free") + theme_bw() +
        labs(title="Demo A: 3-D HSGP (X,Y,Z_bathy) latent-field recovery",
             subtitle=sprintf("vS2 truth; basis %dx%dx%d (Z capped)", MX,MY,MZ),
             x="true log_lambda", y="posterior mean log_lambda"))
# length-scale recovery
print(ggplot(recovery %>% filter(param %in% c("lx","ly","lz")),
             aes(factor(species), mean, ymin=q025, ymax=q975, colour=param)) +
        geom_pointrange(position=position_dodge(0.4)) +
        geom_point(aes(y=true), shape=4, size=3, position=position_dodge(0.4)) +
        theme_bw() + labs(title="Length-scale recovery (x = truth)",
                          x="species", y="length scale (km / m)"))
invisible(dev.off())
cat(sprintf("\nSaved %s and hsgp3D_recovery.pdf\n", file.path(OUTPUT_DIR, "hsgp3D_fit.rds")))
