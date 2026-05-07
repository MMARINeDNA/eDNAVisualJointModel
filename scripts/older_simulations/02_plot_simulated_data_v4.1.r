# =============================================================================
# plot_simulated_data_v4.1.R
#
# Three-page summary of the v4.1 simulation (zero-mean GP, 3 species,
# 500 stations x 6 sample depths). With Option A applied, the simulated
# field is a pure zero-mean GP draw - there's no closed-form expected
# lambda surface, but we can recover a continuous version by KRIGING
# the GP draw from station locations onto a 1 km grid (conditional
# mean: f_grid = K_gs * K_ss^-1 * f_stations). This gives the visual
# of v4's continuous lambda surface while staying faithful to the
# specific GP draw the simulation produced.
#
#   Page 1 - Kriged log(lambda) surface on a 1 km grid for each
#            species at z = 0 m, with station locations overlaid as
#            open circles. Same look as v4's page 1.
#   Page 2 - Per-station log(lambda_true) vs each spatial covariate
#            (X, Y, Z_bathy), one row per species.
#   Page 3 - Water-column sampling-depth effect exp(zsample_pref) on
#            the reversed depth axis, one curve per species.
# =============================================================================

library(tidyverse)
library(patchwork)
library(viridis)

sim <- readRDS("outputs/whale_edna_output_v4.1/whale_edna_sim_v4.1.rds")

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
# Bathymetry contour helper (still relevant for showing shelf geometry).
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

sp_colours <- viridis(n_species, option = "viridis", end = 0.85)

# Per-station truth dataframe (long form by species), surface samples
# only for the X-Y map.
surf_idx <- which(as.numeric(samples[["Z_sample"]]) == 0)

truth_long <- bind_rows(lapply(seq_len(n_species), function(s) {
  lambda_vec <- as.numeric(lambda_true[, s])
  tibble(
    species     = sp_common[s],
    sample_idx  = seq_len(nrow(samples)),
    X           = as.numeric(samples[["X"]]),
    Y           = as.numeric(samples[["Y"]]),
    Z_bathy     = as.numeric(samples[["Z_bathy"]]),
    Z_sample    = as.numeric(samples[["Z_sample"]]),
    lambda_true = lambda_vec,
    log_lambda  = log(lambda_vec)
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

# =============================================================================
# Page 1 - 1 km grid of log(lambda), kriged from station-level GP draws,
#          with stations overlaid as open white circles. Matches v4's
#          page 1 layout.
# =============================================================================
cat("Building 1 km grid (",
    X_km_max + 1, " x ", Y_km_max + 1, " cells)...\n", sep = "")

grid_1km <- expand.grid(
  X = seq(0, X_km_max, by = 1),
  Y = seq(0, Y_km_max, by = 1)
) %>%
  mutate(Z_bathy = bathy_mean_fn(X, Y))

# Anisotropic SE covariance helper. Returns matrix of shape (n1, n2).
aniso_cross_cov <- function(c1, c2, sigma, lx, ly, lz) {
  d2 <- outer(c1[, 1], c2[, 1], FUN = function(a, b) ((a - b) / lx)^2) +
        outer(c1[, 2], c2[, 2], FUN = function(a, b) ((a - b) / ly)^2) +
        outer(c1[, 3], c2[, 3], FUN = function(a, b) ((a - b) / lz)^2)
  sigma * sigma * exp(-0.5 * d2)
}

# Conditioning data: station-level GP coords. Multiple samples at the
# same station share the same (X, Y, Z_bathy), so we collapse to one
# row per station. The simulated field is identical (up to 1e-6 jitter)
# across those redundant rows anyway.
station_first_idx <- match(seq_len(nrow(stations)), samples$station)
station_coords    <- as.matrix(samples[station_first_idx,
                                       c("X", "Y", "Z_bathy")])
gp_field_si       <- sim$truth$gp_field_si        # N samples x n_species
station_f         <- gp_field_si[station_first_idx, , drop = FALSE]
N_st              <- nrow(station_coords)

grid_coords <- as.matrix(grid_1km[, c("X", "Y", "Z_bathy")])
N_grid      <- nrow(grid_coords)

cat(sprintf("Kriging GP draws onto 1 km grid for %d species...\n", n_species))
cat(sprintf("  conditioning on %d stations -> %d grid cells\n", N_st, N_grid))

batch_size <- 25000L

for (s in seq_len(n_species)) {
  p <- gp_params[[s]]
  t0 <- Sys.time()

  # K_ss + jitter so chol() is numerically stable
  K_ss <- aniso_cross_cov(station_coords, station_coords,
                          p$sigma, p$lx, p$ly, p$lz) +
          diag(1e-4, N_st)
  chol_K <- chol(K_ss)
  alpha  <- backsolve(chol_K,
                      forwardsolve(t(chol_K), station_f[, s]))

  f_grid <- numeric(N_grid)
  for (b in seq(1L, N_grid, by = batch_size)) {
    idx  <- b:min(b + batch_size - 1L, N_grid)
    K_gs <- aniso_cross_cov(grid_coords[idx, , drop = FALSE],
                            station_coords,
                            p$sigma, p$lx, p$ly, p$lz)
    f_grid[idx] <- as.numeric(K_gs %*% alpha)
  }
  grid_1km[[paste0("log_lambda_s", s)]] <- p$mu + f_grid
  cat(sprintf("  %s: %.1f s\n",
              sp_common[s], as.numeric(Sys.time() - t0,
                                       units = "secs")))
}

# Shared colour limits across species (1st and 99th percentiles), so
# species panels are visually comparable on the log-lambda scale.
log_lam_all <- unlist(lapply(seq_len(n_species),
                             function(s) grid_1km[[paste0("log_lambda_s", s)]]))
clim <- quantile(log_lam_all, c(0.01, 0.99))

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
      subtitle = expression("Kriged " * log(lambda) * " (z = 0 m, 1 km grid)"),
      x        = "Easting (km, 0 = 100000 UTM E)",
      y        = "Northing (km, 0 = SF, 1270 = 49 deg N)"
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
    title    = "Simulated eDNA fields - v4.1 (zero-mean GP)",
    subtitle = expression("Kriged " * log(lambda) *
                          " on a 1 km grid; open circles = station locations; white contours = 200 m / 1000 m isobaths"),
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9,  colour = "grey40"))
  ) &
  theme(legend.position = "right")

# =============================================================================
# Page 2 - Per-station log(lambda) vs each spatial covariate
# =============================================================================
covariate_levels <- c(
  "X (km, easting)",
  "Y (km, northing)",
  "Z_bathy (m, bottom depth)"
)

page2_df <- bind_rows(
  truth_long %>% transmute(species, covariate = covariate_levels[1],
                           axis_val = X, log_lambda = log_lambda),
  truth_long %>% transmute(species, covariate = covariate_levels[2],
                           axis_val = Y, log_lambda = log_lambda),
  truth_long %>% transmute(species, covariate = covariate_levels[3],
                           axis_val = Z_bathy, log_lambda = log_lambda)
) %>%
  mutate(covariate = factor(covariate, levels = covariate_levels))

page2 <- ggplot(page2_df, aes(x = axis_val, y = log_lambda, colour = species)) +
  geom_point(alpha = 0.4, size = 0.9) +
  facet_grid(species ~ covariate, scales = "free", switch = "y") +
  scale_colour_viridis_d(option = "D", end = 0.85, guide = "none") +
  labs(
    title    = "Per-station truth vs spatial covariates (v4.1)",
    subtitle = "log(lambda) at simulated stations - field is a single zero-mean GP draw per species",
    x        = NULL,
    y        = expression(log(lambda))
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title       = element_text(size = 12, face = "bold"),
        plot.subtitle    = element_text(size = 9,  colour = "grey40"),
        panel.grid.minor = element_blank(),
        strip.placement  = "outside")

# =============================================================================
# Page 3 - Water-column sampling-depth effect exp(zsample_pref) per species
# =============================================================================
zsample_df <- bind_rows(lapply(seq_len(n_species), function(s) {
  p <- gp_params[[s]]
  tibble(
    species    = sp_common[s],
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
    geom_point(colour = col, size = 3.0) +
    geom_text(aes(label = sprintf("%+.2f", log_offset)),
              hjust = -0.30, size = 3.0, colour = "grey25") +
    scale_y_reverse(breaks = sample_depths,
                    labels = paste0(sample_depths, " m")) +
    scale_x_continuous(limits = xlim_mult, expand = c(0.05, 0)) +
    labs(
      title    = sp_common[s],
      x        = expression(exp(delta[Z_sample])),
      y        = if (s == 1) "Sample depth (m)" else NULL
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "italic", size = 11),
          axis.title.y  = if (s == 1) element_text() else element_blank(),
          axis.text.y   = if (s == 1) element_text() else element_blank(),
          axis.ticks.y  = if (s == 1) element_line() else element_blank())
})

page3 <- wrap_plots(p_row3, nrow = 1) +
  plot_annotation(
    title    = "eDNA water column effect - v4.1",
    subtitle = "Per-species multiplicative shift on eDNA at each sample depth",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9,  colour = "grey40"))
  )

# =============================================================================
# Write multi-page PDF
# =============================================================================
out_pdf <- "outputs/whale_edna_output_v4.1/simulated_edna_fields_v4.1.pdf"
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
grDevices::cairo_pdf(out_pdf, width = 11, height = 9, onefile = TRUE)
print(page1)
print(page2)
print(page3)
invisible(dev.off())

cat(sprintf("Saved %s (3 pages)\n", out_pdf))
