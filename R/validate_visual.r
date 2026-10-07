# =============================================================================
# validate_visual.r
#
# Multi-replicate recovery study for stan/hsgp_visual.stan: simulate several
# independent datasets (one field + sightings per seed), fit each species, and
# summarise convergence, field R^2 and 95% CI coverage across replicates.
#
#   Rscript R/validate_visual.r
# Env: SEEDS (e.g. "1:5" or "1,2,3"), CHAINS, WARMUP, SAMPLE, ADAPT_DELTA,
#      TREEDEPTH, OUT (output dir, default outputs/visual/validation)
# Out: <OUT>/validation_{runs,recovery,floor,summary}.csv
# =============================================================================

source("R/fit_visual.r")

env_int <- function(k, d) as.integer(Sys.getenv(k, d))
SEEDS <- eval(parse(text = sprintf("c(%s)", Sys.getenv("SEEDS", "1:5"))))
OUT   <- Sys.getenv("OUT", "outputs/visual/validation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

mod  <- compile_visual()
runs <- list(); recs <- list(); floors <- list()
for (seed in SEEDS) {
  sim <- simulate_visual(seed = seed)
  for (sp in sim$meta$species) {
    r <- fit_visual(sim, sp, mod = mod,
                    chains = env_int("CHAINS", "4"), warmup = env_int("WARMUP", "1000"),
                    sample = env_int("SAMPLE", "1000"),
                    adapt_delta = as.numeric(Sys.getenv("ADAPT_DELTA", "0.9")),
                    max_treedepth = env_int("TREEDEPTH", "12"), seed = 42L + seed)
    runs[[length(runs) + 1]] <- data.frame(
      seed = seed, species = sp, n_det = r$n_det, runtime_min = r$runtime_min,
      divergences = r$divergences, treedepth_hits = r$treedepth_hits,
      min_ebfmi = min(r$ebfmi), max_rhat = r$max_rhat, M = r$M,
      floor_flags = sum(r$floor$flag), R2 = r$field$R2,
      rmse = r$field$rmse, bias = r$field$bias)
    recs[[length(recs) + 1]] <- cbind(seed = seed, r$recovery)
    if (!is.null(r$floor)) floors[[length(floors) + 1]] <- cbind(seed = seed, r$floor)
    cat(sprintf("seed %d %-8s n_det=%4d  %.1f min  div=%d  maxRhat=%.3f  R2=%.3f\n",
                seed, sp, r$n_det, r$runtime_min, r$divergences, r$max_rhat, r$field$R2))
  }
}
runs <- do.call(rbind, runs); recs <- do.call(rbind, recs)
write.csv(runs, file.path(OUT, "validation_runs.csv"), row.names = FALSE)
write.csv(recs, file.path(OUT, "validation_recovery.csv"), row.names = FALSE)
if (length(floors)) write.csv(do.call(rbind, floors), file.path(OUT, "validation_floor.csv"), row.names = FALSE)

cov <- aggregate(covered ~ species + param, data = recs[!is.na(recs$covered), ],
                 FUN = function(x) sprintf("%d/%d", sum(x), length(x)))
rel <- aggregate(cbind(rel_bias = (median - truth) / abs(truth)) ~ species + param,
                 data = recs[!is.na(recs$truth), ], FUN = median)
summary <- merge(cov, rel)
write.csv(summary, file.path(OUT, "validation_summary.csv"), row.names = FALSE)

options(width = 200)
cat("\n=== Per-fit diagnostics ===\n");  print(runs, digits = 3, row.names = FALSE)
cat("\n=== Coverage (truth in 95% CI) and median relative bias across replicates ===\n")
print(summary, digits = 3, row.names = FALSE)
