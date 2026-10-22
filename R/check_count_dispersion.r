# =============================================================================
# check_count_dispersion.r
#
# Do clustered sightings make the Poisson visual model's length-scales follow
# the prior floor, as in the real fits (explainer, "Refit at twice the basis
# budget")? ROADMAP_2026-10-21 section 3 item 1.
#
# Simulates sightings on the REAL 35N effort (data/processed/muri_visual.rds:
# 1,641 segments with their lengths and bottom depths), from a field whose
# length-scales are well above the floor of the real design (LX, LY; the
# floor at the 800-function budget is about 75 km x 106 km), with Poisson or
# negative-binomial groups per segment (NB_SIZES; Inf = Poisson), and fits
# the visual model to each. If clustering alone drives the floor-following,
# the Poisson fits to NB data should flag the floor and the fits to Poisson
# data should not.
#
#   Rscript R/check_count_dispersion.r
# Env: SEEDS ("1"), SPECIES ("humpback,pwsd"), NB_SIZES ("Inf,1"),
#      LX, LY (truth length-scales, km; 150, 250), MU_PWSD (truth log mean
#      density for PWSD; log(0.07), so detections are near the real 89),
#      CHAINS, WARMUP, SAMPLE (4 x 500/500), M_MAX (basis budget; 800),
#      REAL (1: also fit the real 35N data, the reference; default 1),
#      OUT (outputs/visual/count_dispersion)
# Out: <OUT>/fit_<species>_k<size>_s<seed>.rds (summary per fit; existing
#      ones are skipped), <OUT>/count_dispersion_{runs,floor,recovery}.csv
# =============================================================================

source("R/fit_visual.r")

env_int <- function(k, d) as.integer(Sys.getenv(k, d))
env_vec <- function(k, d) trimws(strsplit(Sys.getenv(k, d), ",")[[1]])
SEEDS    <- eval(parse(text = sprintf("c(%s)", Sys.getenv("SEEDS", "1"))))
SPECIES  <- env_vec("SPECIES", "humpback,pwsd")
NB_SIZES <- as.numeric(env_vec("NB_SIZES", "Inf,1"))
LX <- as.numeric(Sys.getenv("LX", "150")); LY <- as.numeric(Sys.getenv("LY", "250"))
MU_PWSD <- as.numeric(Sys.getenv("MU_PWSD", as.character(log(0.07))))
M_MAX <- env_int("M_MAX", as.character(GP_BASIS_BUDGET))
OUT <- Sys.getenv("OUT", "outputs/visual/count_dispersion")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

real_seg <- readRDS("data/processed/muri_visual.rds")$design$segments
gp <- default_gp_params()
for (sp in c("humpback", "pwsd")) { gp[[sp]]$lx <- LX; gp[[sp]]$ly <- LY }
gp$pwsd$mu <- MU_PWSD

mod <- compile_visual()

# Reference: the same model on the real 35N sightings (no truth; floor check only)
if (env_int("REAL", "1") == 1L) {
  vis <- readRDS("data/processed/muri_visual.rds")
  for (sp in SPECIES) {
    f <- file.path(OUT, sprintf("fit_%s_real.rds", sp))
    if (file.exists(f)) next
    sd <- format_stan_data_visual(vis, sp, M_max = M_MAX)
    t0 <- Sys.time()
    fit <- mod$sample(data = sd, chains = env_int("CHAINS", "4"), parallel_chains = env_int("CHAINS", "4"),
                      iter_warmup = env_int("WARMUP", "500"), iter_sampling = env_int("SAMPLE", "500"),
                      adapt_delta = 0.9, max_treedepth = 12L, seed = 42L, init = visual_init(sd),
                      refresh = 0, show_messages = FALSE)
    fl <- cbind(species = sp, floor_check(fit$draws(c("gp_l[1]", "gp_l[2]"), format = "draws_matrix"),
                                          attr(sd, "gp_prior")))
    cnt <- sd$seg_count
    r <- list(setting = data.frame(seed = NA, species = sp, nb_size = NA, lx = NA, ly = NA, n_det = sd$n,
                                   var_mean = var(cnt) / mean(cnt)),
              runtime_min = as.numeric(Sys.time() - t0, units = "mins"),
              divergences = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent),
              max_rhat = max(fit$summary(c("mu_sp", "gp_sigma", "gp_l_raw", "log_sigma"), "rhat")$rhat),
              floor = fl, gp = fit$summary(c("gp_sigma[1]", "gp_l[1]", "gp_l[2]", "D_mean")))
    saveRDS(r, f)
    cat(sprintf("REAL %-8s n_det=%d | %.1f min div=%d maxRhat=%.3f | l = %s (floor %s) flags %d\n", sp, sd$n,
                r$runtime_min, r$divergences, r$max_rhat, paste(round(fl$post_median), collapse = ", "),
                paste(round(fl$ell_min), collapse = ", "), sum(fl$flag)))
  }
}

res <- list()
for (seed in SEEDS) for (k in NB_SIZES) {
  sim <- simulate_visual(seed = seed, gp_params = gp, segments = real_seg,
                         nb_size = setNames(rep(k, 2), c("humpback", "pwsd")))
  for (sp in SPECIES) {
    f <- file.path(OUT, sprintf("fit_%s_k%s_s%d.rds", sp, format(k), seed))
    if (file.exists(f)) { res[[f]] <- readRDS(f); next }
    cnt <- sim$observed[[sp]]$seg_count
    r <- fit_visual(sim, sp, mod = mod, M_max = M_MAX,
                    chains = env_int("CHAINS", "4"), warmup = env_int("WARMUP", "500"),
                    sample = env_int("SAMPLE", "500"), seed = 42L + seed)
    r$fit <- NULL
    r$setting <- data.frame(seed = seed, species = sp, nb_size = k, lx = LX, ly = LY,
                            n_det = r$n_det, var_mean = var(cnt) / mean(cnt))
    saveRDS(r, f); res[[f]] <- r
    cat(sprintf("seed %d %-8s k=%-4s n_det=%4d var/mean %.2f | %.1f min div=%d maxRhat=%.3f R2=%.3f | l = %s (floor %s) flags %d\n",
                seed, sp, format(k), r$n_det, var(cnt) / mean(cnt), r$runtime_min, r$divergences, r$max_rhat,
                r$field$R2, paste(round(r$floor$post_median), collapse = ", "),
                paste(round(r$floor$ell_min), collapse = ", "), sum(r$floor$flag)))
  }
}

real <- lapply(Filter(file.exists, file.path(OUT, sprintf("fit_%s_real.rds", SPECIES))), readRDS)
if (length(real)) {
  cat("\n=== Real 35N data (Poisson model) ===\n")
  print(do.call(rbind, lapply(real, function(r) cbind(r$setting[, c("species", "n_det", "var_mean")], r$floor[, -1]))),
        digits = 3, row.names = FALSE)
}
runs <- do.call(rbind, lapply(res, function(r) cbind(r$setting,
  runtime_min = r$runtime_min, divergences = r$divergences, max_rhat = r$max_rhat,
  R2 = r$field$R2, rmse = r$field$rmse, bias = r$field$bias, floor_flags = sum(r$floor$flag))))
floor <- do.call(rbind, lapply(res, function(r) cbind(r$setting[, c("seed", "species", "nb_size")], r$floor[, -1])))
recov <- do.call(rbind, lapply(res, function(r) cbind(r$setting[, c("seed", "nb_size")], r$recovery)))
write.csv(runs,  file.path(OUT, "count_dispersion_runs.csv"), row.names = FALSE)
write.csv(floor, file.path(OUT, "count_dispersion_floor.csv"), row.names = FALSE)
write.csv(recov, file.path(OUT, "count_dispersion_recovery.csv"), row.names = FALSE)
options(width = 200)
print(runs, digits = 3, row.names = FALSE); print(floor, digits = 3, row.names = FALSE)
