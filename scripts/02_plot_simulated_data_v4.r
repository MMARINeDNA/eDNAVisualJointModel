# =============================================================================
# plot_simulated_data_v4.R
#
# Three-panel figure for the V4 simulation (SF to US/Canada border, rotated
# bathymetry, 300 stations × 3 sample depths with variable per-sample
# replication). Output is a multi-page PDF with one page per panel:
#
#   Page 1 — Expected density surface λ(X, Y, Z_bathy(X, Y)) evaluated on a
#             1 km grid (no kriging — GP mean evaluated in closed form at
#             every grid cell). Station locations overlaid.
#
#   Page 2 — 3 × 3 grid of the true marginal relationship between each
#             species' expected density λ and each of the three spatial
#             covariates (X, Y, Z_bathy). True relationship only (line,
#             no points). Rows = species, columns = covariate; each row
#             has its own y-scale because densities differ by orders of
#             magnitude between species.
#
#   Page 3 — Water column sampling depth effect exp(δ_Z_sample) on the
#             reversed depth axis (0, 150, 500 m).
#
# Lambda is rendered via expression() / plotmath so the Greek letter
# always shows up, and the final PDF is written with cairo_pdf for
# correct Unicode support.
# =============================================================================

library(tidyverse)
library(patchwork)
library(viridis)

sim <- readRDS("outputs/whale_edna_sim_v4.rds")

sp_common      <- sim$meta$sp_common
n_species      <- sim$meta$n_species
sample_depths  <- sim$meta$sample_depths
D              <- sim$meta$n_sample_depth
gp_params      <- sim$truth$gp_params
lambda_true    <- sim$truth$lambda_true_si
samples        <- sim$design$samples
stations       <- sim$design$stations
X_km_max       <- sim$meta$X_km_max
Y_km_max       <- sim$meta$Y_km_max
bathy          <- sim$meta$bathy_profile

# ---------------------------------------------------------------------------
# Re-create the deterministic functions used in the sim (rotated bathymetry
# and Gaussian preferences) from parameters stored in `sim$meta` / gp_params.
# ---------------------------------------------------------------------------
rot <- bathy$rotation_deg * pi / 180

X_prime_fn <- function(X_km, Y_km) {
  raw <- cos(rot) * (X_km / X_km_max) + sin(rot) * (Y_km / Y_km_max)
  300 * raw / (cos(rot) + sin(rot))
}

bathy_mean_fn <- function(X_km, Y_km) {
  xp <- X_prime_fn(X_km, Y_km)
  pmax(
    bathy$abyssal_depth +
      (bathy$shelf_depth - bathy$abyssal_depth) /
      (1 + exp(-bathy$slope_k * (xp - bathy$shelf_break_xprime))),
    10
  )
}

gauss_pref <- function(x, mu, sd, amp) {
  raw <- dnorm(x, mu, sd) / dnorm(mu, mu, sd)
  amp * (raw - 0.5)
}

# Expected log(λ) at (X_km, Y_km) for species s — GP mean only, no draw.
expected_log_lambda <- function(X_km, Y_km, Z_bathy, p) {
  p$mu +
    gauss_pref(Z_bathy, p$zbathy_pref_mu, p$zbathy_pref_sd, p$zbathy_pref_amp) +
    gauss_pref(Y_km,    p$y_pref_mu,      p$y_pref_sd,      p$y_pref_amp)
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
  mutate(Z_bathy = bathy_mean_fn(X, Y))

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
    scale_x_continuous(breaks = seq(0, X_km_max, 200)) +
    scale_y_continuous(breaks = seq(0, Y_km_max, 200)) +
    labs(
      title    = sp_common[s],
      subtitle = expression("Expected " * lambda * " (z = 0 m, 1 km grid)"),
      x        = "Easting (km, 0 = 100000 UTM E)",
      y        = "Northing (km, 0 = SF, 1270 = 49°N)"
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
    title    = "Simulated eDNA fields - SF to US/Canada border (v4)",
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
# For each species we hold the "other" variables at representative values:
#
#   λ vs X     :  Y fixed at the species' preferred latitude y_pref_mu;
#                 Z_bathy is computed from (X, Y_fixed) via the rotated
#                 shelf-slope profile, so the X response folds in the
#                 bathymetric preference.
#
#   λ vs Y     :  X fixed at the domain centre (X_km_max / 2); Z_bathy
#                 is computed from (X_fixed, Y). Combines the latitude
#                 preference and the Z_bathy-via-Y gradient.
#
#   λ vs Z_bathy : Y fixed at y_pref_mu; X irrelevant since the mean
#                 function uses Z_bathy directly. This isolates the
#                 bathymetric preference curve.
#
# Free y-scales per species (densities differ by orders of magnitude)
# and free x-scales per covariate (X, Y, Z_bathy have different ranges).
# =============================================================================

covariate_levels <- c(
  "X (km, easting)",
  "Y (km, northing)",
  "Z_bathy (m, bottom depth)"
)
cov_X <- covariate_levels[1]
cov_Y <- covariate_levels[2]
cov_Z <- covariate_levels[3]

# Build the sweeps for all species × covariates
row2_df <- bind_rows(lapply(seq_len(n_species), function(s) {
  p <- gp_params[[s]]

  # X sweep: Y fixed at species' preferred latitude
  x_seq    <- seq(0, X_km_max, length.out = 300)
  Y_fixed  <- p$y_pref_mu
  zb_x     <- bathy_mean_fn(x_seq, Y_fixed)
  ll_x     <- expected_log_lambda(x_seq, Y_fixed, zb_x, p)

  # Y sweep: X fixed at domain centre
  y_seq    <- seq(0, Y_km_max, length.out = 300)
  X_fixed  <- X_km_max / 2
  zb_y     <- bathy_mean_fn(X_fixed, y_seq)
  ll_y     <- expected_log_lambda(X_fixed, y_seq, zb_y, p)

  # Z_bathy sweep: Y fixed at species' preferred latitude
  z_seq    <- seq(10, 2600, length.out = 300)
  ll_z     <- expected_log_lambda(NA_real_, Y_fixed, z_seq, p)

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

page2 <- ggplot(row2_df, aes(x = axis_val, y = lambda, colour = species)) +
  geom_line(linewidth = 1.1) +
  facet_grid(species ~ covariate, scales = "free") +
  scale_colour_manual(values = sp_colours, guide = "none") +
  labs(
    title    = "True marginal density relationships - v4",
    subtitle = expression(
      "Expected " * lambda * " (animals/km"^2 *
      ") swept over one covariate at a time; others held at representative values"
    ),
    x        = NULL,
    y        = expression(lambda ~ "(animals/km"^2 * ")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    plot.subtitle   = element_text(size = 9,  colour = "grey40"),
    strip.text.x    = element_text(size = 10),
    strip.text.y    = element_text(size = 10, face = "italic"),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# Page 3 — Water column sampling depth effect (unchanged from v3/v4)
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
    title    = "eDNA water column effect - v4",
    subtitle = "Per-species multiplicative shift on eDNA at each sample depth",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  )

# =============================================================================
# Write multi-page PDF (one page per panel) with cairo_pdf for Unicode
# =============================================================================

out_pdf <- "outputs/simulated_edna_fields_v4.pdf"
grDevices::cairo_pdf(out_pdf, width = 8, height = 8, onefile = TRUE)
print(page1)
print(page2)
print(page3)
invisible(dev.off())

cat(sprintf("Saved %s (3 pages)\n", out_pdf))
