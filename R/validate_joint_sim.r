# =============================================================================
# validate_joint_sim.r
#
# Roadmap Phase 4 check: the joint simulator is consistent with BOTH
# single-source models. For each seed, simulate_joint() draws one field; the
# eDNA part is fitted with HSGP4eDNA's gp2d_bathysp model (R/fit_edna.r) and
# the visual part with stan/hsgp_visual.stan (R/fit_visual.r). Each must still
# recover the shared field.
#
#   Rscript R/validate_joint_sim.r
# Env: SEEDS (default "1:2"), OUT (default outputs/joint_sim/validation)
# Out: <OUT>/joint_sim_{edna,visual}.csv
# =============================================================================

source("R/functions_joint.R")
source("R/fit_edna.r")
source("R/fit_visual.r")

SEEDS <- eval(parse(text = sprintf("c(%s)", Sys.getenv("SEEDS", "1:2"))))
OUT   <- Sys.getenv("OUT", "outputs/joint_sim/validation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

vmod <- compile_visual()
edna_rows <- list(); vis_rows <- list(); agree_rows <- list()
for (seed in SEEDS) {
  sim <- simulate_joint(seed = seed)
  agree_rows[[length(agree_rows) + 1]] <- cbind(seed = seed, joint_field_agreement(sim))

  # eDNA half: HSGP4eDNA's validated bathysp_surface settings, matched spline df
  scen <- sprintf("bathysp_joint_s%d", seed)
  e <- fit_edna(scen, sim = sim$edna, basis = c(17L, 6L), df = 4L,
                chains = 3L, warmup = 400L, sample = 400L,
                log_file = file.path(OUT, sprintf("edna_s%d.log", seed)))
  edna_rows[[length(edna_rows) + 1]] <- data.frame(
    seed = seed, species = e$field$species, R2 = e$field$R2, rmse = e$field$rmse,
    divergences = sum(e$diagnostic_summary$num_divergent),
    max_rhat = max(e$recovery$rhat, na.rm = TRUE), runtime_min = e$runtime_min)
  cat(sprintf("seed %d eDNA: R2 %s | div %d | max Rhat %.3f | %.0f min\n", seed,
              paste(sprintf("%.2f", e$field$R2), collapse = " / "),
              sum(e$diagnostic_summary$num_divergent), max(e$recovery$rhat, na.rm = TRUE),
              e$runtime_min))

  # Visual half
  for (sp in sim$visual$meta$species) {
    r <- fit_visual(sim$visual, sp, mod = vmod, seed = 42L + seed)
    rec <- r$recovery
    vis_rows[[length(vis_rows) + 1]] <- data.frame(
      seed = seed, species = sp, n_det = r$n_det, R2 = r$field$R2,
      divergences = r$divergences, max_rhat = r$max_rhat, floor_flags = sum(r$floor$flag),
      D_covered = rec$covered[rec$param == "D_mean"],
      sigma_covered = rec$covered[rec$param == "sigma_c"],
      gp_covered = sprintf("%d/3", sum(rec$covered[rec$param %in% c("gp_sigma[1]", "gp_l[1]", "gp_l[2]")])),
      runtime_min = r$runtime_min)
    cat(sprintf("seed %d visual %-8s n_det %4d: R2 %.2f | div %d | max Rhat %.3f\n",
                seed, sp, r$n_det, r$field$R2, r$divergences, r$max_rhat))
  }
}
edna <- do.call(rbind, edna_rows); vis <- do.call(rbind, vis_rows); agree <- do.call(rbind, agree_rows)
write.csv(edna, file.path(OUT, "joint_sim_edna.csv"), row.names = FALSE)
write.csv(vis,  file.path(OUT, "joint_sim_visual.csv"), row.names = FALSE)
write.csv(agree, file.path(OUT, "joint_sim_field_agreement.csv"), row.names = FALSE)
options(width = 200)
cat("\n=== Shared field: station/segment pairs < 10 km ===\n"); print(agree, digits = 3, row.names = FALSE)
cat("\n=== eDNA half (HSGP4eDNA gp2d_bathysp) ===\n"); print(edna, digits = 3, row.names = FALSE)
cat("\n=== Visual half (hsgp_visual) ===\n");          print(vis, digits = 3, row.names = FALSE)
