# =============================================================================
# tests/smoke.R
#
# Fast end-to-end check that every model in this repo still compiles, samples,
# and recovers its simulated truth. Run from the repo root before committing
# any change under R/, stan/, or a submodule bump:
#
#   Rscript tests/smoke.R
#
# Small data, small basis, short chains: the thresholds are deliberately loose.
# A failure means "something is broken", not "the model is mis-calibrated";
# the full validation runs (README results table) answer the latter.
#
# Cases:
#   edna    HSGP4eDNA gp2d_bathysp (2-D HSGP + bottom-depth spline) via
#           R/fit_edna.r, on a small bathysp_surface-style simulation
#   visual  stan/hsgp_visual.stan (2-D HSGP + bottom-depth spline, line
#           transects) via R/fit_visual.r, both species of simulate_visual()
#
# Logs and outputs go to outputs/smoke/ (git-ignored). Exit status is non-zero
# if any check fails.
# =============================================================================

suppressMessages({
  library(readr)
})

SMOKE_DIR <- "outputs/smoke"
dir.create(SMOKE_DIR, showWarnings = FALSE, recursive = TRUE)

results <- list()
check <- function(case, what, ok, value) {
  results[[length(results) + 1L]] <<- data.frame(
    case = case, check = what, value = value, pass = isTRUE(ok))
}
timed <- function(expr) {
  t0 <- Sys.time(); force(expr)
  as.numeric(Sys.time() - t0, units = "mins")
}

# -----------------------------------------------------------------------------
# Case 1: eDNA (HSGP4eDNA gp2d_bathysp)
# -----------------------------------------------------------------------------
cat("=== [edna] HSGP4eDNA gp2d_bathysp ===\n")
source("R/fit_edna.r")
edna_log <- file.path(SMOKE_DIR, "edna.log")
edna <- NULL
mins <- timed(edna <- tryCatch(
  fit_edna("bathysp_smoke",
           sim_args = list(seed = 1L, n_stations = 100L,
                           sample_depths = 0, zsample_pref = NULL),
           basis = c(8L, 4L), df = 4L,
           chains = 2L, warmup = 200L, sample = 200L,
           log_file = edna_log),
  error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL }))
cat(sprintf("  %.1f min (log: %s)\n", mins, edna_log))

check("edna", "ran", !is.null(edna), if (is.null(edna)) "error" else "ok")
if (!is.null(edna)) {
  ndiv <- sum(edna$diagnostic_summary$num_divergent)
  rhat <- max(edna$recovery$rhat, na.rm = TRUE)
  r2   <- setNames(edna$field$R2, edna$field$species)
  check("edna", "divergences == 0", ndiv == 0, ndiv)
  check("edna", "max Rhat (GP hypers) < 1.1", rhat < 1.1, round(rhat, 3))
  check("edna", "hake field R2 > 0.5", r2[["Pacific hake"]] > 0.5,
        round(r2[["Pacific hake"]], 2))
}

# -----------------------------------------------------------------------------
# Case 2: visual (hsgp_visual: 2-D HSGP + bottom-depth spline, line transects)
# -----------------------------------------------------------------------------
cat("=== [visual] hsgp_visual ===\n")
source("R/fit_visual.r")
vis <- list()
mins <- timed({
  vsim <- simulate_visual(seed = 1L)
  vmod <- tryCatch(compile_visual(quiet = TRUE),
                   error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL })
  if (!is.null(vmod)) for (sp in vsim$meta$species) {
    vis[[sp]] <- tryCatch(
      fit_visual(vsim, sp, mod = vmod, chains = 2L, warmup = 300L, sample = 300L,
                 HSGP_M = c(8L, 4L)),
      error = function(e) { message("  ERROR (", sp, "): ", conditionMessage(e)); NULL })
  }
})
cat(sprintf("  %.1f min\n", mins))

for (sp in c("humpback", "pwsd")) {
  r <- vis[[sp]]
  check("visual", sprintf("%s ran", sp), !is.null(r), if (is.null(r)) "error" else "ok")
  if (is.null(r)) next
  check("visual", sprintf("%s divergences == 0", sp), r$divergences == 0, r$divergences)
  check("visual", sprintf("%s max Rhat < 1.1", sp), r$max_rhat < 1.1, round(r$max_rhat, 3))
}
if (!is.null(vis$humpback)) {
  r  <- vis$humpback
  sc <- r$recovery[r$recovery$param == "sigma_c", ]
  check("visual", "humpback field R2 > 0.5", r$field$R2 > 0.5, round(r$field$R2, 2))
  check("visual", "humpback sigma truth in 95% CI", isTRUE(sc$covered),
        sprintf("%.3g in [%.3g, %.3g]", sc$truth, sc$q025, sc$q975))
}

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
res <- do.call(rbind, results)
cat("\n=== Smoke test results ===\n")
options(width = 200)
print(res, row.names = FALSE, right = FALSE)
write_csv(res, file.path(SMOKE_DIR, "smoke_results.csv"))
if (!all(res$pass)) {
  cat(sprintf("\nSMOKE TEST FAILED: %d of %d checks\n", sum(!res$pass), nrow(res)))
  quit(status = 1)
}
cat(sprintf("\nSMOKE TEST PASSED (%d checks)\n", nrow(res)))
