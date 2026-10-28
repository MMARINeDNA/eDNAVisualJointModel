# =============================================================================
# fit_joint_real.r
#
# The joint model on the REAL data (ROADMAP_2026-10-21, section 5 item 7):
# the eDNA and visual objects from R/01_load_muri.r (data/processed/, tracked;
# no private files needed), put together by real_joint_object(), and fitted as
#   edna        eDNA only: the real eDNA fit (the joint model with the visual
#               likelihood off equals HSGP4eDNA's own eDNA model, and shares
#               the joint fit's basis and priors);
#   visual      visual only;
#   joint_conv  both sources, the cetaceans' conversion factors estimated
#               (the real conversion factors are unknown: the eDNA object's
#               are simulator placeholders);
#   joint       both sources with the placeholder conversion factors fixed
#               (not meaningful for the cetaceans on real data; for checks).
# eDNA options: Ct with per-plate curves estimated from the standards, the
# dilution offsets and per-bottle volumes (from the object), the sampling-depth
# smooth over the target depths. Visual options come from the object
# (truth$sp_params): humpback hazard-rate + group size; PWSD half-normal +
# group size + Beaufort + visibility; negative-binomial counts for both.
#
#   MODES=edna Rscript R/fit_joint_real.r
# Env: MODES (default "edna,visual,joint_conv"), QPCR_CURVE ("estimated" or
#      "fixed"), ZDEPTH ("smooth" or "none"), CT_SD ("expected" or
#      "conditional"), JUNK_NOISE_DIST (2 non-centred, 1 centred), CONV_FREE
#      (default "humpback,pwsd"), M_MAX (basis budget; 800), CHAINS, WARMUP,
#      SAMPLE (4 x 1000/1000), ADAPT_DELTA (0.9), TREEDEPTH (12), REFRESH
#      (100), OUT (outputs/joint/real), CSV_DIR (keep the CmdStan chain CSVs
#      there; default: a temporary folder, deleted at the end), GRID_KM
#      (prediction grid over water, data/processed/pred_grid_real_<km>km.rds
#      from R/make_pred_grid.r; default 20, 0 = none), PPC (1 = predictive-
#      check summaries, default)
# Out: <OUT>/fit_<mode>.rds (summaries: field at stations and segments with
#      90% intervals, parameters, each chain's means, floor check,
#      diagnostics, the field on the grid, predictive checks), and
#      <OUT>/real_{params,floor,runs}.csv over every mode found; with a grid,
#      <OUT>/real_maps{,_sd}.png; with PPC, ppc_{visual,edna}_<mode>.png.
# =============================================================================

source("R/fit_joint.r")

env_int <- function(k, d) as.integer(Sys.getenv(k, d))
MODES <- strsplit(Sys.getenv("MODES", "edna,visual,joint_conv"), ",")[[1]]
CONV_FREE <- strsplit(Sys.getenv("CONV_FREE", "humpback,pwsd"), ",")[[1]]
OPTS <- list(qpcr_curve = Sys.getenv("QPCR_CURVE", "estimated"), zdepth = Sys.getenv("ZDEPTH", "smooth"),
             ct_sd = Sys.getenv("CT_SD", "expected"), junk_noise_dist = env_int("JUNK_NOISE_DIST", "2"),
             M_max = env_int("M_MAX", as.character(GP_BASIS_BUDGET)))
OUT <- Sys.getenv("OUT", "outputs/joint/real")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CSV_DIR <- Sys.getenv("CSV_DIR", "")
GRID_KM <- as.numeric(Sys.getenv("GRID_KM", "20"))
PPC <- env_int("PPC", "1")
grid <- if (GRID_KM > 0) pred_grid_real(GRID_KM)
if (nzchar(CSV_DIR)) dir.create(CSV_DIR, showWarnings = FALSE, recursive = TRUE)

real <- real_joint_object()
cat(sprintf("real data: %d bottles (%d stations), %d qPCR wells, %d standards, %d MB replicates; %d segments\n",
            nrow(real$edna$design$samples), length(unique(real$edna$design$samples$station)),
            length(real$edna$observed$qpcr_sample_idx), length(real$edna$observed$qpcr_std_plate_idx),
            nrow(real$edna$observed$mb_reads), nrow(real$visual$design$segments)))

if (length(MODES) && any(nzchar(MODES))) {
  mod <- compile_joint()
  for (m in MODES) {
    cat(sprintf("\n=== real %s (qPCR curves %s, depth %s, junk noise %s) ===\n", m, OPTS$qpcr_curve,
                OPTS$zdepth, c("centred", "non-centred")[OPTS$junk_noise_dist]))
    r <- do.call(fit_joint, c(list(real, if (m == "joint_conv") "joint" else m,
                                   conv_free = if (m == "joint_conv") CONV_FREE else character(0),
                                   mod = mod, chains = env_int("CHAINS", "4"),
                                   warmup = env_int("WARMUP", "1000"), sample = env_int("SAMPLE", "1000"),
                                   adapt_delta = as.numeric(Sys.getenv("ADAPT_DELTA", "0.9")),
                                   max_treedepth = env_int("TREEDEPTH", "12"),
                                   refresh = env_int("REFRESH", "100"), seed = 42L,
                                   output_dir = if (nzchar(CSV_DIR)) CSV_DIR,
                                   pred = grid, use_ppc = PPC), OPTS))
    r$mode <- m; r$opts <- OPTS; r$pred <- grid
    r$params$mode <- m; r$floor$mode <- m          # fit_joint() labels joint_conv "joint"
    cat(sprintf("%s: %.0f min | divergences %d | treedepth %d | max Rhat %.3f | E-BFMI %s\n",
                m, r$runtime_min, r$divergences, r$treedepth_hits, r$max_rhat,
                paste(sprintf("%.2f", r$ebfmi), collapse = ", ")))
    slow <- r$params$variable[which(r$params$rhat > 1.05)]
    if (length(slow)) { cat("chain means where R-hat > 1.05:\n"); print(r$chain_means[c("lp__", slow), ], digits = 3) }
    r$fit <- NULL; r$stan_data <- NULL
    saveRDS(r, file.path(OUT, sprintf("fit_%s.rds", m)))
  }
}

# ---- Tables over every mode found ----------------------------------------------
files <- list.files(OUT, "^fit_.*\\.rds$", full.names = TRUE)
files <- files[basename(files) != "fit_prior.rds"]              # the prior predictive: not in the tables
if (!length(files)) quit(save = "no")
res <- lapply(files, readRDS)
runs <- do.call(rbind, lapply(res, function(r) data.frame(mode = r$mode, runtime_min = r$runtime_min,
  divergences = r$divergences, treedepth_hits = r$treedepth_hits, min_ebfmi = min(r$ebfmi), max_rhat = r$max_rhat)))
params <- do.call(rbind, lapply(res, `[[`, "params"))
floor  <- do.call(rbind, lapply(res, `[[`, "floor"))
write.csv(runs,   file.path(OUT, "real_runs.csv"), row.names = FALSE)
write.csv(params, file.path(OUT, "real_params.csv"), row.names = FALSE)
write.csv(floor,  file.path(OUT, "real_floor.csv"), row.names = FALSE)
options(width = 200)
cat("\n=== Runs ===\n"); print(runs, digits = 3, row.names = FALSE)
cat("\n=== Main parameters ===\n")
print(params[!grepl("^(beta_bathy|alpha_ct_plate|beta_ct_plate|beta_zdepth)", params$variable), ],
      digits = 3, row.names = FALSE)
cat("\n=== Length-scale floor check ===\n"); print(floor, digits = 3, row.names = FALSE)
cat(sprintf("\nSaved %s/fit_<mode>.rds and real_{runs,params,floor}.csv\n", OUT))

# ---- Maps and predictive checks ----------------------------------------------------
source("R/plot_joint.r")
names(res) <- vapply(res, `[[`, "", "mode")
gres <- Filter(function(r) !is.null(r$grid), res[intersect(c("edna", "visual", "joint_conv", "joint"), names(res))])
if (length(gres)) {
  pg <- gres[[1]]$pred
  ht <- 3.4 * length(gres[[1]]$grid) + 0.8
  ggsave(file.path(OUT, "real_maps.png"),
         plot_field_maps(gres, pg, real, modes = names(gres),
                         title = sprintf("Real data: posterior mean log density on a %g km grid", GRID_KM),
                         subtitle = "Dots: eDNA stations; orange: transect segments (Beaufort <= 5)"),
         width = 2.1 * length(gres) + 1.2, height = ht, dpi = 160, bg = "white")
  ggsave(file.path(OUT, "real_maps_sd.png"),
         plot_field_maps(gres, pg, real, what = "sd", modes = names(gres),
                         title = "Real data: posterior SD of log density"),
         width = 2.1 * length(gres) + 1.2, height = ht, dpi = 160, bg = "white")
  cat(sprintf("Saved %s/real_maps{,_sd}.png\n", OUT))
}
if (file.exists(file.path(OUT, "fit_prior.rds"))) {
  post <- res[[intersect(c("joint_conv", "joint", "visual", "edna"), names(res))[1]]]
  ggsave(file.path(OUT, "prior_pred.png"),
         plot_prior_pred(readRDS(file.path(OUT, "fit_prior.rds")), post,
                         sprintf("Real data, prior predictive; blue: %s posterior predictive", MODE_LABEL[[post$mode]])),
         width = 10, height = 8, dpi = 150, bg = "white")
}
for (m in names(res)) {
  r <- res[[m]]
  if (!is.null(r$ppc$visual))
    ggsave(file.path(OUT, sprintf("ppc_visual_%s.png", m)),
           plot_ppc_visual(r, sprintf("Real data, posterior predictive checks, visual: %s fit", MODE_LABEL[[m]])),
           width = 13, height = 3.3 * length(r$ppc$visual) + 0.8, dpi = 150, bg = "white")
  if (!is.null(r$ppc$edna))
    ggsave(file.path(OUT, sprintf("ppc_edna_%s.png", m)),
           plot_ppc_edna(r, sprintf("Real data, posterior predictive checks, eDNA: %s fit", MODE_LABEL[[m]])),
           width = 10, height = 6.5, dpi = 150, bg = "white")
}
