# =============================================================================
# plot_simulated_data_v2.R
#
# Three-panel figure for the V2 simulation (Oregon / Washington coast,
# bathymetry as a function of X only, 150 stations × 5 sample depths).
# Output is a multi-page PDF with one page per panel:
#
#   Page 1 — Expected density surface λ on a 1 km grid (closed-form GP
#             mean evaluated at every grid cell). Station locations
#             overlaid.
#
#   Page 2 — 3 × 3 grid of the true marginal relationship between each
#             species' expected density λ and each spatial covariate
#             (X, Y, Z_bathy). True relationship only — no points, no
#             second geom. Rows = species, columns = covariate; each
#             species column has its own y-scale.
#
#             V2 has no latitude preference and bathymetry depends only
#             on X, so the λ vs Y panels are flat by construction —
#             they're kept in the layout for parity with later versions.
#
#   Page 3 — Water column sampling depth effect exp(δ_Z_sample) on the
#             reversed depth axis (V2 uses 5 depths: 0, 50, 150, 300,
#             500 m).
#
# Lambda is rendered via expression() / plotmath so the Greek letter
# always shows up, and the final PDF is written with cairo_pdf for
# correct Unicode support.
# =============================================================================

library(tidyverse)
library(patchwork)
library(viridis)

sim <- readRDS("outputs/whale_edna_output_v2/whale_edna_sim_v2.rds")

sp_common      <- sim$meta$sp_common
n_species      <- sim$meta$n_species
sample_depths  <- sim$meta$sample_depths
D              <- sim$meta$n_sample_depth
gp_params      <- sim$truth$gp_params
lambda_true    <- sim$truth$lambda_true_si
samples        <- sim$design$samples
stations       <- sim$design$stations

# v2 doesn't store these in meta — derive them from the UTM box.
X_km_max <- (sim$meta$X_max_utm - sim$meta$X_min_utm) / 1000   # 300 km
Y_km_max <- (sim$meta$Y_max_utm - sim$meta$Y_min_utm) / 1000   # 400 km

# ---------------------------------------------------------------------------
# Re-create the deterministic functions used in the v2 sim. Bathymetry is
# a function of X alone (no rotation, no Y dependence). There is no
# latitude preference (no y_pref* fields in gp_params).
# ---------------------------------------------------------------------------
bathy_mean_fn <- function(X_km, Y_km = NULL) {
  abyssal <- 2500
  shelf   <- 80
  pmax(
    abyssal + (shelf - abyssal) / (1 + exp(-0.06 * (X_km - 180))),
    10
  )
}

gauss_pref <- function(x, mu, sd, amp) {
  raw <- dnorm(x, mu, sd) / dnorm(mu, mu, sd)
  amp * (raw - 0.5)
}

# Expected log(λ) at (X, Y, Z_bathy) for species s — GP mean only.
# v2's mean function is mu + zbathy_pref(Z_bathy); no Y dependence.
expected_log_lambda <- function(X_km, Y_km, Z_bathy, p) {
  p$mu +
    gauss_pref(Z_bathy, p$zbathy_pref_mu, p$zbathy_pref_sd, p$zbathy_pref_amp)
}

sp_colours <- viridis(n_species, option = "viridis", end = 0.85)

# =============================================================================
# Page 1 — Expected density surface, 1 km grid
# =============================================================================

cat("Building 1 km grid (",
    X_km_max + 1, " × ", Y_km_max + 1, " cells)...\n", sep = "")

grid_1km <- expand.grid(
  X = seq(0, X_km_max, by = 1),
  Y = seq(0, Y_km_max, by = 1)
) %>%
  mutate(Z_bathy = bathy_mean_fn(X))

for (s in seq_len(n_species)) {
  p <- gp_params[[s]]
  grid_1km[[paste0("log_lambda_s", s)]] <-
    expected_log_lambda(grid_1km$X, grid_1km$Y, grid_1km$Z_bathy, p)
}

# Shared colour limits across species (1st and 99th percentiles on log scale)
log_lam_all <- unlist(lapply(seq_len(n_species),
                             function(s) grid_1km[[paste0("log_lambda_s", s)]]))
clim <- quantile(log_lam_all, c(0.01, 0.99))

# 200 m and ~1000 m isobaths
iso_levels <- c(200, 1000)

p_row1 <- lapply(seq_len(n_species), function(s) {

  grid_df <- grid_1km %>%
    transmute(X, Y, Z_bathy,
              log_lambda = .data[[paste0("log_lambda_s", s)]])

  ggplot() +
    geom_raster(data = grid_df, aes(X, Y, fill = log_lambda)) +
    geom_contour(
      data      = grid_df,
      aes(X, Y, z = Z_bathy),
      breaks    = iso_levels,
      colour    = "white",
      linewidth = 0.25,
      alpha     = 0.5
    ) +
    geom_point(
      data   = stations,
      aes(X, Y),
      shape  = 21,
      size   = 0.8,
      stroke = 0.3,
      fill   = NA,
      colour = "white"
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = clim,
      name   = expression(log(lambda)),
      oob    = scales::squish
    ) +
    coord_fixed(
      ratio  = 1,
      xlim   = c(0, X_km_max),
      ylim   = c(0, Y_km_max),
      expand = FALSE
    ) +
    scale_x_continuous(breaks = seq(0, X_km_max, 100)) +
    scale_y_continuous(breaks = seq(0, Y_km_max, 100)) +
    labs(
      title    = sp_common[s],
      subtitle = expression("Expected " * lambda * " (z = 0 m, 1 km grid)"),
      x        = "Easting (km from offshore boundary)",
      y        = "Northing (km from Cape Blanco)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title      = element_text(face = "italic", size = 11),
      plot.subtitle   = element_text(size = 8.5, colour = "grey40"),
      legend.position = if (s == n_species) "right" else "none",
      axis.title.y    = if (s == 1) element_text() else element_blank(),
      axis.text.y     = if (s == 1) element_text() else element_blank(),
      axis.ticks.y    = if (s == 1) element_line() else element_blank(),
      panel.grid      = element_blank()
    )
})

page1 <- wrap_plots(p_row1, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Simulated eDNA fields - Oregon / Washington coast (v2)",
    subtitle = expression("Expected density " * lambda *
                          " on a 1 km grid at z = 0 m"),
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  ) &
  theme(legend.position = "right")

# =============================================================================
# Page 2 — True marginal relationships (3 species × 3 covariates)
#
# v2 mean function is mu + zbathy_pref(Z_bathy); Z_bathy depends only on
# X. So:
#
#   λ vs X     :  Z_bathy = bathy(X), expected λ = exp(mu + zbathy_pref(.)).
#                 The X response folds in the bathymetric preference.
#   λ vs Y     :  λ does not depend on Y in v2; line is flat by design.
#                 Kept in the grid for parity with v3+ plots.
#   λ vs Z_bathy : exp(mu + zbathy_pref(Z_bathy)). Standalone curve.
#
# Free y-scales per species (densities differ); free x-scales per covariate.
# =============================================================================

covariate_levels <- c(
  "X (km, easting)",
  "Y (km, northing)",
  "Z_bathy (m, bottom depth)"
)
cov_X <- covariate_levels[1]
cov_Y <- covariate_levels[2]
cov_Z <- covariate_levels[3]

row2_df <- bind_rows(lapply(seq_len(n_species), function(s) {
  p <- gp_params[[s]]

  # X sweep — Z_bathy varies through bathy_mean_fn(X)
  x_seq    <- seq(0, X_km_max, length.out = 300)
  zb_x     <- bathy_mean_fn(x_seq)
  ll_x     <- expected_log_lambda(x_seq, NA_real_, zb_x, p)

  # Y sweep — λ doesn't depend on Y. Hold X at domain centre and Z_bathy
  # at the corresponding shelf-slope value, so the line shows the
  # constant exp(mu + zbathy_pref(Z_bathy(X_centre))).
  y_seq    <- seq(0, Y_km_max, length.out = 300)
  X_fixed  <- X_km_max / 2
  zb_y     <- bathy_mean_fn(rep(X_fixed, length(y_seq)))
  ll_y     <- expected_log_lambda(X_fixed, y_seq, zb_y, p)

  # Z_bathy sweep — direct
  z_seq    <- seq(10, 2600, length.out = 300)
  ll_z     <- expected_log_lambda(NA_real_, NA_real_, z_seq, p)

  bind_rows(
    data.frame(species = sp_common[s], covariate = cov_X,
               axis_val = x_seq, lambda = exp(ll_x)),
    data.frame(species = sp_common[s], covariate = cov_Y,
               axis_val = y_seq, lambda = exp(ll_y)),
    data.frame(species = sp_common[s], covariate = cov_Z,
               axis_val = z_seq, lambda = exp(ll_z))
  )
})) %>%
  mutate(
    species   = factor(species,   levels = sp_common),
    covariate = factor(covariate, levels = covariate_levels)
  )

# Species in columns, covariates in rows, with per-species y-scales.
page2_species <- lapply(seq_len(n_species), function(s) {
  df <- row2_df %>% filter(species == sp_common[s])
  p <- ggplot(df, aes(x = axis_val, y = lambda)) +
    geom_line(colour = sp_colours[s], linewidth = 1.1) +
    facet_wrap(~ covariate, nrow = 3, scales = "free_x",
               strip.position = "right") +
    labs(title = sp_common[s], x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "italic", size = 11, hjust = 0.5),
      strip.text.y     = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
  if (s == 1) {
    p <- p + labs(y = expression(lambda ~ "(animals/km"^2 * ")"))
  }
  p
})

page2 <- wrap_plots(page2_species, nrow = 1) +
  plot_annotation(
    title    = "True marginal density relationships - v2",
    subtitle = expression(
      "Expected " * lambda * " (animals/km"^2 *
      "); v2 has no Y dependence so middle row is flat by design"
    ),
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  )

# =============================================================================
# Page 3 — Water column sampling depth effect (5 depths in v2)
# =============================================================================

zsample_df <- bind_rows(lapply(seq_len(n_species), function(s) {
  p <- gp_params[[s]]
  tibble(
    species    = sp_common[s],
    depth_idx  = seq_len(D),
    Z_sample   = sample_depths,
    log_offset = p$zsample_pref,
    edna_mult  = exp(p$zsample_pref)
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

xlim_mult <- c(
  min(0.05, min(zsample_df$edna_mult) * 0.85),
  max(2.00, max(zsample_df$edna_mult) * 1.12)
)

p_row3 <- lapply(seq_len(n_species), function(s) {
  df  <- filter(zsample_df, species == sp_common[s])
  col <- sp_colours[s]

  ggplot(df, aes(x = edna_mult, y = Z_sample)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    geom_path(colour = col, linewidth = 1.0) +
    geom_point(colour = col, size = 3.5) +
    geom_text(
      aes(label = sprintf("%+.2f", log_offset)),
      hjust = -0.30, size = 3.2, colour = "grey25"
    ) +
    scale_y_reverse(breaks = sample_depths,
                    labels = paste0(sample_depths, " m")) +
    scale_x_continuous(limits = xlim_mult, expand = c(0.05, 0)) +
    labs(
      title    = sp_common[s],
      subtitle = "eDNA water column effect",
      x        = expression(exp(delta[Z_sample])),
      y        = "Sample depth (m)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "italic", size = 11),
      plot.subtitle = element_text(size = 8.5, colour = "grey40"),
      axis.title.y  = if (s == 1) element_text() else element_blank(),
      axis.text.y   = if (s == 1) element_text() else element_blank(),
      axis.ticks.y  = if (s == 1) element_line() else element_blank()
    )
})

page3 <- wrap_plots(p_row3, nrow = 1) +
  plot_annotation(
    title    = "eDNA water column effect - v2",
    subtitle = "Per-species multiplicative shift on eDNA at each sample depth",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  )

# =============================================================================
# Write multi-page PDF (one page per panel) with cairo_pdf for Unicode
# =============================================================================

out_pdf <- "outputs/whale_edna_output_v2/simulated_edna_fields_v2.pdf"
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
grDevices::cairo_pdf(out_pdf, width = 8, height = 8, onefile = TRUE)
print(page1)
print(page2)
print(page3)
invisible(dev.off())

cat(sprintf("Saved %s (3 pages)\n", out_pdf))
