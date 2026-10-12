# Simulate one visual survey (seed 1, 25 transects) and fit both species with the spatial
# visual model. Run from the repo root: Rscript docs/explainer/02_fit_visual_sim.R
source("docs/explainer/00_common.R")
source("R/fit_visual.r")
sim <- simulate_visual(seed = 1L)          # simulate before loading anything else
saveRDS(sim, file.path(EXPL_OUT, "visual_sim.rds"))
mod <- compile_visual()
for (sp in c("humpback", "pwsd")) {
  sd <- format_stan_data_visual(sim, sp)
  cat(sprintf("%s: n=%d n_seg=%d M=%d\n", sp, sd$n, sd$n_seg, sd$M))
  saveRDS(sd, file.path(EXPL_OUT, sprintf("visual_sim_%s_standata.rds", sp)))
  t0 <- Sys.time()
  fit <- mod$sample(data = sd, chains = 3, parallel_chains = 3, iter_warmup = 800, iter_sampling = 800,
                    adapt_delta = 0.9, max_treedepth = 12, seed = 42L, init = visual_init(sd), refresh = 100)
  cat("runtime (min):", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "\n")
  fit$save_object(file.path(EXPL_OUT, sprintf("visual_sim_%s_fit.rds", sp)))
  print(fit$diagnostic_summary())
}
