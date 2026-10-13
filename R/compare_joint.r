# =============================================================================
# compare_joint.r
#
# Roadmap Phase 5 check: on one simulate_joint() dataset, fit
# stan/hsgp_joint.stan as eDNA-only, visual-only and joint (same model, priors
# and basis; only the data differ), plus a joint fit with the cetaceans' eDNA
# conversion factors estimated. Compares field recovery (R^2, 95% CI width and
# coverage of the true log density at every station and segment) and
# hyperparameter recovery.
#
#   Rscript R/compare_joint.r
# Env: SEED (default 1), MODES (default "edna,visual,joint,joint_conv"),
#      CONV_FREE (species freed in joint_conv; default "humpback,pwsd"),
#      CHAINS, WARMUP, SAMPLE (default 4 x 800/800), OUT (default
#      outputs/joint/compare_s<SEED>)
# Each mode's summary is saved as <OUT>/fit_<mode>.rds, so modes can run in
# separate processes; the tables and figure are (re)built from whatever
# summaries exist.
# Out: <OUT>/compare_{field,recovery,runs}.csv, <OUT>/compare_joint.{pdf,png}
# =============================================================================

source("R/fit_joint.r")
suppressMessages(library(ggplot2))

env_int <- function(k, d) as.integer(Sys.getenv(k, d))
SEED  <- env_int("SEED", "1")
MODES <- strsplit(Sys.getenv("MODES", "edna,visual,joint,joint_conv"), ",")[[1]]
CONV_FREE <- strsplit(Sys.getenv("CONV_FREE", "humpback,pwsd"), ",")[[1]]
OUT <- Sys.getenv("OUT", sprintf("outputs/joint/compare_s%d", SEED))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
MODE_ORDER <- c("edna", "visual", "joint", "joint_conv")

sim <- simulate_joint(seed = SEED)
# modes of one seed may run as parallel processes: write via a rename so a
# reader never sees a half-written file
tmp <- tempfile("sim", tmpdir = OUT, fileext = ".rds")
saveRDS(sim, tmp); file.rename(tmp, file.path(OUT, "sim.rds"))
if (length(MODES) && any(nzchar(MODES))) {
  mod <- compile_joint()
  for (m in MODES) {
    cat(sprintf("\n=== %s ===\n", m))
    r <- fit_joint(sim, if (m == "joint_conv") "joint" else m,
                   conv_free = if (m == "joint_conv") CONV_FREE else character(0),
                   mod = mod, chains = env_int("CHAINS", "4"),
                   warmup = env_int("WARMUP", "800"), sample = env_int("SAMPLE", "800"),
                   seed = 42L + SEED)
    r$field$mode <- m; r$recovery$mode <- m; r$mode <- m
    cat(sprintf("%s: %.0f min | divergences %d | treedepth %d | max Rhat %.3f | E-BFMI %s\n",
                m, r$runtime_min, r$divergences, r$treedepth_hits, r$max_rhat,
                paste(sprintf("%.2f", r$ebfmi), collapse = ", ")))
    print(r$field, digits = 3, row.names = FALSE)
    r$fit <- NULL; r$stan_data <- NULL
    saveRDS(r, file.path(OUT, sprintf("fit_%s.rds", m)))
  }
}

# ---- Assemble whatever modes have finished ------------------------------------
files <- file.path(OUT, sprintf("fit_%s.rds", MODE_ORDER))
res <- setNames(lapply(files[file.exists(files)], readRDS), MODE_ORDER[file.exists(files)])
if (!length(res)) quit(save = "no")
field <- do.call(rbind, lapply(res, `[[`, "field"))
# A species with no data in a mode (hake in visual-only) has a prior-only field
lt <- sim$meta$lt_species
field <- field[!(field$mode == "visual" & !field$species %in% lt), ]
recov <- do.call(rbind, lapply(res, `[[`, "recovery"))
runs  <- do.call(rbind, lapply(res, function(r) data.frame(
  mode = r$mode, runtime_min = r$runtime_min, divergences = r$divergences,
  treedepth_hits = r$treedepth_hits, min_ebfmi = min(r$ebfmi), max_rhat = r$max_rhat)))
write.csv(field, file.path(OUT, "compare_field.csv"), row.names = FALSE)
write.csv(recov, file.path(OUT, "compare_recovery.csv"), row.names = FALSE)
write.csv(runs,  file.path(OUT, "compare_runs.csv"), row.names = FALSE)
options(width = 200)
cat("\n=== Runs ===\n");  print(runs, digits = 3, row.names = FALSE)
cat("\n=== Field recovery (all stations + segments) ===\n")
print(field[order(field$species, match(field$mode, MODE_ORDER)), ], digits = 3, row.names = FALSE)
cat("\n=== Hyperparameter / conversion-factor recovery ===\n")
print(recov[order(recov$variable, match(recov$mode, MODE_ORDER)),
            c("mode", "variable", "truth", "median", "q025", "q975", "covered", "rhat")],
      digits = 3, row.names = FALSE)

# ---- Figure: truth vs posterior mean, per species x mode ----------------------
truth <- sim$field$log_lambda_si; src <- sim$field$locations$source
sp <- colnames(res[[1]]$post_mean)
pts <- do.call(rbind, lapply(names(res), function(m) do.call(rbind, lapply(seq_along(sp), function(s)
  data.frame(mode = m, species = sp[s], source = src, truth = truth[, s],
             post = res[[m]]$post_mean[, s])))))
pts <- pts[!(pts$mode == "visual" & !pts$species %in% lt), ]
pts$mode <- factor(pts$mode, MODE_ORDER[MODE_ORDER %in% names(res)])
lab <- field; lab$label <- sprintf("R2 %.2f\nCI %.2f", lab$R2, lab$ci_width)
lab$mode <- factor(lab$mode, levels(pts$mode))
g <- ggplot(pts, aes(truth, post, colour = source)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.4) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.1, vjust = 1.2, size = 3) +
  facet_grid(species ~ mode, scales = "free") +
  geom_text(data = data.frame(mode = factor("visual", levels(pts$mode)),
                              species = setdiff(sp, lt), label = "no visual data"),
            aes(x = 0, y = 0, label = label), inherit.aes = FALSE, colour = "grey40", size = 3.5) +
  scale_colour_manual(values = c(station = "#D55E00", segment = "#0072B2")) +
  labs(x = "true log density", y = "posterior mean log density", colour = NULL,
       title = sprintf("Joint vs single-source field recovery (simulate_joint seed %d)", SEED),
       subtitle = "R2 and mean 95% CI width of log density over all stations and segments (hake has no visual data)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "compare_joint.pdf"), g, width = 10, height = 7.5)
ggsave(file.path(OUT, "compare_joint.png"), g, width = 10, height = 7.5, dpi = 150)
cat(sprintf("\nSaved %s/compare_{field,recovery,runs}.csv and compare_joint.{pdf,png}\n", OUT))
