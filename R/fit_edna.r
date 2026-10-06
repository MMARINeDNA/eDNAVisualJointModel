# =============================================================================
# fit_edna.r
#
# Thin driver for the eDNA-only HSGP model. The model, simulator, and Stan-data
# formatter all live in the HSGP4eDNA submodule (external/HSGP4eDNA); nothing
# is copied here. This file only runs them:
#
#   1. simulate_bathysp() (HSGP4eDNA R/functions.R) -> outputs/<scenario>/sim.rds
#   2. R/09_fit_gp2d_bathysp.r  (2-D HSGP + bottom-depth spline,
#      stan/hsgp_2d_bathysp.stan) -> outputs/<scenario>/fit_gp2d_bathysp_*
#
# Both steps run with the submodule root as working directory, because the
# HSGP4eDNA scripts use repo-relative paths. Outputs therefore land in
# external/HSGP4eDNA/outputs/<scenario>/ (git-ignored by that repo).
#
# Usage (from this repo's root):
#   source("R/fit_edna.r")
#   res <- fit_edna("bathysp_surface")             # full validation settings
#   res$field                                      # latent-field R^2 per species
# =============================================================================

HSGP4EDNA_DIR <- "external/HSGP4eDNA"

# Run `expr` with the HSGP4eDNA submodule root as working directory.
in_hsgp4edna <- function(expr) {
  if (!file.exists(file.path(HSGP4EDNA_DIR, "R", "functions.R")))
    stop("HSGP4eDNA submodule missing: run `git submodule update --init`")
  old <- setwd(HSGP4EDNA_DIR)
  on.exit(setwd(old))
  force(expr)
}

# scenario : name of the outputs/ subdirectory. Must start with "bathysp" so
#            HSGP4eDNA's 09 driver treats the fit as matched (spline truth).
# sim_args : passed to simulate_bathysp(); defaults reproduce
#            HSGP4eDNA's 01_sim_bathysp_surface.r.
# basis    : HSGP basis per axis (MX, MY); df: spline df for Z_bathy.
# log_file : where the fit's console output goes (NULL = this console).
fit_edna <- function(scenario  = "bathysp_surface",
                     sim_args  = list(seed = 202L, sample_depths = 0,
                                      zsample_pref = NULL),
                     basis     = c(17L, 6L),
                     df        = 5L,
                     chains    = 3L,
                     warmup    = 400L,
                     sample    = 400L,
                     log_file  = NULL) {
  stopifnot(grepl("^bathysp", scenario), length(basis) == 2L)
  # Resolve before changing directory into the submodule.
  # (normalizePath() leaves a not-yet-existing file relative, so resolve the
  # directory and re-attach the file name.)
  if (!is.null(log_file)) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    log_file <- file.path(normalizePath(dirname(log_file)), basename(log_file))
  }

  in_hsgp4edna({
    out_dir <- file.path("outputs", scenario)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # 1. Simulate in a clean environment so HSGP4eDNA's helpers do not leak
    #    into the caller's session.
    env <- new.env()
    sys.source("R/functions.R", envir = env)
    sim <- do.call(env$simulate_bathysp, sim_args)
    saveRDS(sim, file.path(out_dir, "sim.rds"))

    # 2. Fit in a subprocess, exactly as HSGP4eDNA runs it.
    status <- system2(
      "Rscript", "R/09_fit_gp2d_bathysp.r",
      env    = c(sprintf("SCENARIO=%s", scenario),
                 sprintf("MX=%d", basis[1]), sprintf("MY=%d", basis[2]),
                 sprintf("DF=%d", df), sprintf("CHAINS=%d", chains),
                 sprintf("WARMUP=%d", warmup), sprintf("SAMPLE=%d", sample)),
      stdout = if (is.null(log_file)) "" else log_file,
      stderr = if (is.null(log_file)) "" else log_file
    )
    if (status != 0) stop("HSGP4eDNA 09_fit_gp2d_bathysp.r failed (exit ", status, ")")

    res <- readRDS(file.path(out_dir, "fit_gp2d_bathysp_summary.rds"))
    res$out_dir <- file.path(HSGP4EDNA_DIR, out_dir)
    res
  })
}
