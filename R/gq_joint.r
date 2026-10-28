# =============================================================================
# gq_joint.r
#
# Adds the prediction grid and the predictive checks to a finished real-data
# fit whose CmdStan chain CSVs were kept (fit_joint_real.r with CSV_DIR, e.g.
# the jreal_*_ctcond jobs), without refitting: CmdStan's standalone generated
# quantities run the current model on the saved draws, with the grid and
# use_ppc = 1 in the data. The fit's summary (<OUT>/fit_<mode>.rds) gains
# $grid, $ppc and $pred, and the maps and checks are drawn as
# fit_joint_real.r does.
#
#   MODE=edna CSV_DIR=outputs/joint/real_ctcond/csv_edna OUT=outputs/joint/real_ctcond \
#     Rscript R/gq_joint.r
# Env: MODE (edna, visual, joint, joint_conv), CSV_DIR, OUT, GRID_KM (default
#      20), CONV_FREE (default "humpback,pwsd"); the model options must match
#      the fit's (read from its saved $opts).
# The fits made before the hazard-rate bound used b unbounded, so the data
# here set b_min = 0 (the bound changes nothing else in the model).
# =============================================================================
source("R/fit_joint.r")
source("R/plot_joint.r")

MODE <- Sys.getenv("MODE"); CSV_DIR <- Sys.getenv("CSV_DIR"); OUT <- Sys.getenv("OUT")
stopifnot(MODE %in% c("edna", "visual", "joint", "joint_conv"), dir.exists(CSV_DIR))
GRID_KM <- as.numeric(Sys.getenv("GRID_KM", "20"))
CONV_FREE <- strsplit(Sys.getenv("CONV_FREE", "humpback,pwsd"), ",")[[1]]
rds <- file.path(OUT, sprintf("fit_%s.rds", MODE))
r <- readRDS(rds)
opts <- r$opts
csvs <- list.files(CSV_DIR, "\\.csv$", full.names = TRUE)
csvs <- csvs[!grepl("-(profile|diagnostic)", csvs)]
cat(sprintf("%s: %d chain CSVs in %s; options %s\n", MODE, length(csvs), CSV_DIR,
            paste(names(opts), unlist(opts), sep = "=", collapse = ", ")))

real <- real_joint_object()
grid <- pred_grid_real(GRID_KM)
sw <- MODE_SWITCHES[[if (MODE == "joint_conv") "joint" else MODE]]
sd <- format_stan_data_joint(real, use_edna = sw[1], use_visual = sw[2],
                             conv_free = if (MODE == "joint_conv") CONV_FREE else character(0),
                             qpcr_curve = opts$qpcr_curve, zdepth = opts$zdepth, ct_sd = opts$ct_sd,
                             junk_noise_dist = opts$junk_noise_dist, M_max = opts$M_max,
                             pred = grid, use_ppc = 1L)
sd$b_min <- as.array(rep(0, sd$L))

mod <- compile_joint()
t0 <- Sys.time()
gq <- mod$generate_quantities(fitted_params = csvs, data = sd, parallel_chains = length(csvs))
cat(sprintf("generated quantities: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))

# One draws object with the fit's parameters and transformed parameters and
# the new generated quantities, behind the $draws() interface the summaries use
D <- posterior::bind_draws(cmdstanr::read_cmdstan_csv(csvs, format = "draws_array")$post_warmup_draws,
                           posterior::subset_draws(gq$draws(format = "draws_array"),
                                                   variable = c("log_lambda_pred", "det_pdf", "det_size")),
                           along = "variable")
fit <- list(draws = function(variables = NULL, format = "draws_array") {
  x <- posterior::subset_draws(D, variable = variables)
  switch(format, draws_matrix = posterior::as_draws_matrix(x), x)
})

r$grid <- summarise_grid(fit, sd)
r$ppc  <- ppc_joint(fit, real, sd)
r$pred <- grid
saveRDS(r, rds)
cat(sprintf("Saved the grid and predictive checks into %s\n", rds))

ggsave(file.path(OUT, sprintf("real_maps_%s.png", MODE)),
       plot_field_maps(setNames(list(r), MODE), grid, real, modes = MODE,
                       title = sprintf("Real data, %s fit", MODE_LABEL[[MODE]]),
                       subtitle = "posterior mean log density"),
       width = 4.6, height = 3.4 * sd$S + 0.8, dpi = 160, bg = "white")
if (!is.null(r$ppc$visual))
  ggsave(file.path(OUT, sprintf("ppc_visual_%s.png", MODE)),
         plot_ppc_visual(r, sprintf("Real data, posterior predictive checks, visual: %s fit", MODE_LABEL[[MODE]])),
         width = 13, height = 3.3 * length(r$ppc$visual) + 0.8, dpi = 150, bg = "white")
if (!is.null(r$ppc$edna))
  ggsave(file.path(OUT, sprintf("ppc_edna_%s.png", MODE)),
         plot_ppc_edna(r, sprintf("Real data, posterior predictive checks, eDNA: %s fit", MODE_LABEL[[MODE]])),
         width = 10, height = 6.5, dpi = 150, bg = "white")
cat("Saved the figures in", OUT, "\n")
