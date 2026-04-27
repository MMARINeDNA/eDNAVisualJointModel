# =============================================================================
# plot_simulated_data_v3.2.R
#
# Two-page figure for the v3.2 hake-only, surface-only, zero-mean-GP
# simulation. With Option A applied, the simulated field is a pure
# zero-mean GP draw (no deterministic habitat preference function), so
# there is nothing closed-form to plot on a 1 km grid - the only
# meaningful "truth" lives at the simulated station locations.
#
#   Page 1 - Stations on the X-Y plane, coloured by log(lambda_true).
#            White contours mark the 200 m and 1000 m isobaths so the
#            shelf-slope geometry is still visible.
#
#   Page 2 - log(lambda_true) at each station as a function of the three
#            spatial covariates (X, Y, Z_bathy). Free x-scales per
#            covariate. There's no underlying functional truth to draw
#            a smooth line for - the field is just a single GP draw -
#            so this page shows only the per-station scatter.
# =============================================================================

library(tidyverse)
library(patchwork)
library(viridis)

sim <- readRDS("outputs/whale_edna_output_v3.2/whale_edna_sim_v3.2.rds")

sp_common      <- sim$meta$sp_common
n_species      <- sim$meta$n_species
lambda_true    <- sim$truth$lambda_true_si
samples        <- sim$design$samples
stations       <- sim$design$stations
X_km_max       <- sim$meta$X_km_max
Y_km_max       <- sim$meta$Y_km_max
bathy          <- sim$meta$bathy_profile

stopifnot(n_species == 1L)

# ---------------------------------------------------------------------------
# Re-create the deterministic bathymetry function (still relevant for the
# isobath contours on the map) - the GP mean function is gone with Option A.
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

sp_colour <- viridis(1, option = "viridis", end = 0.85)

# Per-station truth dataframe
lambda_vec <- as.numeric(lambda_true[, 1])
truth_df <- tibble(
  X           = as.numeric(samples[["X"]]),
  Y           = as.numeric(samples[["Y"]]),
  Z_bathy     = as.numeric(samples[["Z_bathy"]]),
  lambda_true = lambda_vec,
  log_lambda  = log(lambda_vec)
)

# =============================================================================
# Page 1 - Stations coloured by simulated log(lambda), with isobaths
# =============================================================================

cat("Building bathymetry contour grid...\n")
grid_iso <- expand.grid(
  X = seq(0, X_km_max, by = 5),
  Y = seq(0, Y_km_max, by = 5)
) %>%
  mutate(Z_bathy = bathy_mean_fn(X, Y))

iso_levels <- c(200, 1000)

page1 <- ggplot() +
  geom_contour(
    data      = grid_iso,
    aes(X, Y, z = Z_bathy),
    breaks    = iso_levels,
    colour    = "grey60",
    linewidth = 0.3,
    alpha     = 0.8
  ) +
  geom_point(
    data   = truth_df,
    aes(X, Y, colour = log_lambda),
    size   = 1.8,
    alpha  = 0.9
  ) +
  scale_colour_viridis_c(
    option = "viridis",
    name   = expression(log(lambda))
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
    title    = sprintf("Simulated eDNA field (v3.2, zero-mean GP) - %s", sp_common[1]),
    subtitle = "Per-station log(lambda) - GP draw with sigma=1.2, (lx, ly, lz)=(50, 300, 150). Grey isobaths: 200, 1000 m.",
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
# Page 2 - log(lambda) at each station vs each spatial covariate
# =============================================================================

covariate_levels <- c(
  "X (km, easting)",
  "Y (km, northing)",
  "Z_bathy (m, bottom depth)"
)

page2_df <- bind_rows(
  truth_df %>% transmute(covariate = covariate_levels[1], axis_val = X,
                         log_lambda = log_lambda),
  truth_df %>% transmute(covariate = covariate_levels[2], axis_val = Y,
                         log_lambda = log_lambda),
  truth_df %>% transmute(covariate = covariate_levels[3], axis_val = Z_bathy,
                         log_lambda = log_lambda)
) %>%
  mutate(covariate = factor(covariate, levels = covariate_levels))

page2 <- ggplot(page2_df, aes(x = axis_val, y = log_lambda)) +
  geom_point(colour = sp_colour, alpha = 0.6, size = 1.2) +
  facet_wrap(~ covariate, ncol = 3, scales = "free_x") +
  labs(
    title    = sprintf("Per-station truth vs spatial covariates (v3.2) - %s", sp_common[1]),
    subtitle = "log(lambda) at simulated stations. No closed-form truth - field is a single GP draw.",
    x        = NULL,
    y        = expression(log(lambda))
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
