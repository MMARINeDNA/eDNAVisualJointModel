# =============================================================================
# plot_simulated_data_V3.R
#
# Three-row figure for the V3 simulation (SF to US/Canada border, rotated
# bathymetry, 300 stations × 3 sample depths):
#
#   Row 1 — True expected animal density λ(X, Y, Z_bathy(X, Y)) evaluated
#            on a 1 km grid (no kriging — the GP mean is evaluated in closed
#            form at every grid cell). Station locations overlaid.
#
#   Row 2 — Marginal habitat response: λ vs Z_bathy, points from stations
#            coloured by latitude Y, line showing the expected curve at
#            each species' preferred Y.
#
#   Row 3 — Water column sampling depth effect exp(δ_Z_sample) on the
#            reversed depth axis (0, 150, 500 m).
#
# =============================================================================

library(tidyverse)
library(patchwork)
library(viridis)

sim <- readRDS("outputs/whale_edna_sim_V3.rds")

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

# Expected log(λ) at (X_km, Y_km) for species s — the GP mean only (no
# random realisation). This is what Row 1 draws on the grid.
expected_log_lambda <- function(X_km, Y_km, Z_bathy, p) {
  p$mu +
    gauss_pref(Z_bathy, p$zbathy_pref_mu, p$zbathy_pref_sd, p$zbathy_pref_amp) +
    gauss_pref(Y_km,    p$y_pref_mu,      p$y_pref_sd,      p$y_pref_amp)
}

sp_colours <- viridis(n_species, option = "viridis", end = 0.85)

# =============================================================================
# Row 1 — Expected density surface, 1 km grid
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

# 200 m and shelf-break (~1290 m here: midpoint of the logistic) isobaths
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
      name   = "log(λ)",
      oob    = scales::squish
    ) +
    coord_fixed(
      ratio = 1,
      xlim  = c(0, X_km_max),
      ylim  = c(0, Y_km_max),
      expand = FALSE
    ) +
    scale_x_continuous(breaks = seq(0, X_km_max, 200)) +
    scale_y_continuous(breaks = seq(0, Y_km_max, 200)) +
    labs(
      title    = sp_common[s],
      subtitle = "Expected λ (z = 0 m, 1 km grid)",
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

# =============================================================================
# Row 2 — λ vs Z_bathy marginal, points coloured by Y
# =============================================================================

# Station-level mean lambda (average across sample depths) for each species
station_lam <- bind_rows(lapply(seq_len(n_species), function(s) {
  data.frame(
    station = samples$station,
    Y       = samples$Y,
    Z_bathy = samples$Z_bathy,
    lambda  = lambda_true[, s],
    species = sp_common[s]
  )
})) %>%
  group_by(station, Y, Z_bathy, species) %>%
  summarise(lambda = mean(lambda), .groups = "drop") %>%
  mutate(species = factor(species, levels = sp_common))

# Expected λ vs Z_bathy at each species' preferred latitude Y = y_pref_mu
zbathy_seq <- seq(10, 2600, by = 5)
bathy_mean_df <- bind_rows(lapply(seq_len(n_species), function(s) {
  p <- gp_params[[s]]
  z_mu <- gauss_pref(zbathy_seq, p$zbathy_pref_mu, p$zbathy_pref_sd, p$zbathy_pref_amp)
  y_mu <- gauss_pref(p$y_pref_mu, p$y_pref_mu,     p$y_pref_sd,      p$y_pref_amp)
  data.frame(
    Z_bathy = zbathy_seq,
    lambda  = exp(p$mu + z_mu + y_mu),
    species = sp_common[s]
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

p_row2 <- ggplot() +
  geom_point(
    data  = station_lam,
    aes(x = Z_bathy, y = lambda, colour = Y),
    size  = 1.3,
    alpha = 0.8
  ) +
  geom_line(
    data      = bathy_mean_df,
    aes(x = Z_bathy, y = lambda),
    colour    = "grey25",
    linewidth = 0.8
  ) +
  facet_wrap(~species, nrow = 1, scales = "free_y") +
  scale_colour_viridis_c(
    option = "mako",
    name   = "Y (km)\nfrom SF",
    limits = c(0, Y_km_max)
  ) +
  scale_x_continuous(breaks = c(0, 200, 500, 1000, 1500, 2000, 2500)) +
  labs(
    title    = "Bathymetric × latitude habitat",
    subtitle = paste("Points = station means (all habitat structure), ",
                     "line = expected λ at the species' preferred Y"),
    x = "Bottom depth Z_bathy (m)",
    y = "λ (animals / km²)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 11),
    plot.subtitle   = element_text(size = 9, colour = "grey40"),
    strip.text      = element_text(face = "italic", size = 10),
    legend.position = "right"
  )

# =============================================================================
# Row 3 — Water column sampling depth effect (3 depths)
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

# =============================================================================
# Assemble
# =============================================================================

row1 <- wrap_plots(p_row1, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

row2 <- p_row2

row3 <- wrap_plots(p_row3, nrow = 1)

fig <- row1 / row2 / row3 +
  plot_layout(heights = c(2.2, 0.9, 0.9)) +
  plot_annotation(
    title    = "Simulated eDNA fields — SF to US/Canada border (V3)",
    subtitle = paste0(
      "Row 1: expected density λ on 1 km grid at z = 0 m; ",
      "Row 2: λ vs Z_bathy (points coloured by Y); ",
      "Row 3: eDNA water column effect at 3 sample depths"
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  )

ggsave("outputs/simulated_edna_fields_V3.png", fig,
       width = 13, height = 18, dpi = 150, limitsize = FALSE)

cat("Saved outputs/simulated_edna_fields_V3.png\n")
