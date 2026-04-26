# =============================================================================
# plot_simulated_data_v3.2.R
#
# Two-page figure for the v3.2 hake-only, surface-only simulation.
#
#   Page 1 - Expected density surface lambda(X, Y, Z_bathy(X, Y)) for
#             hake on a 1 km grid (no kriging - GP mean evaluated in
#             closed form at every grid cell). Station locations
#             overlaid.
#
#   Page 2 - True marginal relationship between hake's expected density
#             lambda and each of the three spatial covariates
#             (X, Y, Z_bathy). True relationship only (line, no points).
#
# The water column / metabarcoding pages from v3 are dropped: surface
# only means the water-column multiplier is 1 everywhere, and there's
# no MB process to plot.
# =============================================================================

library(tidyverse)
library(patchwork)
library(viridis)

sim <- readRDS("outputs/whale_edna_output_v3.2/whale_edna_sim_v3.2.rds")

sp_common      <- sim$meta$sp_common
n_species      <- sim$meta$n_species
gp_params      <- sim$truth$gp_params
lambda_true    <- sim$truth$lambda_true_si
samples        <- sim$design$samples
stations       <- sim$design$stations
X_km_max       <- sim$meta$X_km_max
Y_km_max       <- sim$meta$Y_km_max
bathy          <- sim$meta$bathy_profile

stopifnot(n_species == 1L)

# ---------------------------------------------------------------------------
# Re-create the deterministic functions used in the sim
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

expected_log_lambda <- function(X_km, Y_km, Z_bathy, p) {
  p$mu +
    gauss_pref(Z_bathy, p$zbathy_pref_mu, p$zbathy_pref_sd, p$zbathy_pref_amp) +
    gauss_pref(Y_km,    p$y_pref_mu,      p$y_pref_sd,      p$y_pref_amp)
}

sp_colour <- viridis(1, option = "viridis", end = 0.85)

# =============================================================================
# Page 1 - Expected density surface, 1 km grid
# =============================================================================

cat("Building 1 km grid (",
    X_km_max + 1, " x ", Y_km_max + 1, " cells)...\n", sep = "")

grid_1km <- expand.grid(
  X = seq(0, X_km_max, by = 1),
  Y = seq(0, Y_km_max, by = 1)
) %>%
  mutate(Z_bathy = bathy_mean_fn(X, Y))

p <- gp_params[[1]]
grid_1km$log_lambda <- expected_log_lambda(grid_1km$X, grid_1km$Y, grid_1km$Z_bathy, p)

clim <- quantile(grid_1km$log_lambda, c(0.01, 0.99))

iso_levels <- c(200, 1000)

page1 <- ggplot() +
  geom_raster(data = grid_1km, aes(X, Y, fill = log_lambda)) +
  geom_contour(
    data      = grid_1km,
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
    title    = sprintf("Simulated eDNA field (v3.2) - %s", sp_common[1]),
    subtitle = expression("Expected " * lambda * " on a 1 km grid at z = 0 m (surface)"),
    x        = "Easting (km, 0 = 100000 UTM E)",
    y        = "Northing (km, 0 = SF, 1270 = 49 deg N)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    panel.grid    = element_blank()
  )

# =============================================================================
# Page 2 - True marginal relationships (single species, three covariates)
# =============================================================================

covariate_levels <- c(
  "X (km, easting)",
  "Y (km, northing)",
  "Z_bathy (m, bottom depth)"
)
cov_X <- covariate_levels[1]
cov_Y <- covariate_levels[2]
cov_Z <- covariate_levels[3]

p <- gp_params[[1]]

x_seq    <- seq(0, X_km_max, length.out = 300)
Y_fixed  <- p$y_pref_mu
zb_x     <- bathy_mean_fn(x_seq, Y_fixed)
ll_x     <- expected_log_lambda(x_seq, Y_fixed, zb_x, p)

y_seq    <- seq(0, Y_km_max, length.out = 300)
X_fixed  <- X_km_max / 2
zb_y     <- bathy_mean_fn(X_fixed, y_seq)
ll_y     <- expected_log_lambda(X_fixed, y_seq, zb_y, p)

z_seq    <- seq(10, 2600, length.out = 300)
ll_z     <- expected_log_lambda(NA_real_, Y_fixed, z_seq, p)

row2_df <- bind_rows(
  data.frame(species = sp_common[1], covariate = cov_X,
             axis_val = x_seq, lambda = exp(ll_x)),
  data.frame(species = sp_common[1], covariate = cov_Y,
             axis_val = y_seq, lambda = exp(ll_y)),
  data.frame(species = sp_common[1], covariate = cov_Z,
             axis_val = z_seq, lambda = exp(ll_z))
) %>%
  mutate(covariate = factor(covariate, levels = covariate_levels))

page2 <- ggplot(row2_df, aes(x = axis_val, y = lambda)) +
  geom_line(colour = sp_colour, linewidth = 1.1) +
  facet_wrap(~ covariate, ncol = 3, scales = "free_x") +
  labs(
    title    = sprintf("True marginal density relationships (v3.2) - %s", sp_common[1]),
    subtitle = expression("Expected " * lambda * " (animals/km"^2 *
                          ") swept over one covariate at a time; others held at representative values"),
    x        = NULL,
    y        = expression(lambda ~ "(animals/km"^2 * ")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# Write multi-page PDF
# =============================================================================

out_pdf <- "outputs/whale_edna_output_v3.2/simulated_edna_fields_v3.2.pdf"
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
grDevices::cairo_pdf(out_pdf, width = 8, height = 8, onefile = TRUE)
print(page1)
print(page2)
invisible(dev.off())

cat(sprintf("Saved %s (2 pages)\n", out_pdf))
