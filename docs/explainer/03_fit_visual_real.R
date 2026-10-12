# Spatial visual model on the real MURI line-transect data (both species).
# Run from the repo root: SP=humpback Rscript docs/explainer/03_fit_visual_real.R
# Env: SP (species), M_MAX (basis budget; default GP_BASIS_BUDGET, 800 from HSGP4eDNA v1.4),
#      TAG (output tag; default "real").
# The explainer's first fits (budget 400) were flagged by the floor check and refitted with
# M_MAX=800 TAG=real800; on v1.4 the default is 800, so a doubled refit is M_MAX=1600.
source("docs/explainer/00_common.R")
source("R/fit_visual.r")
sp  <- Sys.getenv("SP", "humpback")
tag <- Sys.getenv("TAG", "real")
M_MAX_ENV <- as.integer(Sys.getenv("M_MAX", as.character(GP_BASIS_BUDGET)))
vis <- readRDS("data/processed/muri_visual.rds")
sd  <- format_stan_data_visual(vis, sp, M_max = M_MAX_ENV)
cat(sprintf("%s: n=%d n_seg=%d w=%.2f det_key=%d use_bft=%d n_bft=%d M=%d K_bathy=%d\n", sp, sd$n,
            sd$n_seg, sd$w, sd$det_key, sd$use_bft_covar, sd$n_bft, sd$M, sd$K_bathy))
saveRDS(sd, file.path(EXPL_OUT, sprintf("visual_%s_%s_standata.rds", tag, sp)))
mod <- compile_visual()
t0 <- Sys.time()
fit <- mod$sample(data = sd, chains = 2, parallel_chains = 2, iter_warmup = 500, iter_sampling = 500,
                  adapt_delta = 0.9, max_treedepth = 12, seed = 42L, init = visual_init(sd), refresh = 50)
cat("runtime (min):", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "\n")
fit$save_object(file.path(EXPL_OUT, sprintf("visual_%s_%s_fit.rds", tag, sp)))
print(fit$diagnostic_summary())
