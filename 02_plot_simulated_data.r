# =============================================================================
# plot_simulated_data.R
#
# Three-row figure:
#   Row 1 — True animal density lambda at the surface, plotted on the X/Y
#            grid at 1 km resolution using the deterministic mean function
#            (mu + zbathy_pref(Z_bathy(X))); sampled locations overlaid.
#
#   Row 2 — Relationship between bathymetry Z_bathy and true animal density
#            lambda. Z_bathy on x-axis, lambda on y-axis; one line per
#            species showing the mean and spread across all stations.
#            This shows the ecological habitat preference over the
#            shelf-slope bathymetric gradient.
#
#   Row 3 — Water column sampling depth effect: the multiplicative effect
#            of Z_sample on eDNA concentration (zsample_effect), shown on
#            a reversed depth axis (shallow at top). This is the separate
#            observation process, distinct from animal density.
#
# =============================================================================

library(tidyverse)
library(patchwork)
library(viridis)
library(dplyr)

sim <- readRDS("whale_edna_sim.rds")

sp_common      <- sim$meta$sp_common
sp_names       <- sim$meta$sp_names
n_species      <- sim$meta$n_species
sample_depths  <- sim$meta$sample_depths    # 0, 50, 150, 300, 500
D              <- sim$meta$n_sample_depth   # 5
gp_params      <- sim$truth$gp_params
lambda_true    <- sim$truth$lambda_true     # N × S  (animal density)
zsample_effect <- sim$truth$zsample_effect  # N × S  (water column factor)
samples        <- sim$design$samples        # X, Y, Z_bathy, Z_sample, etc.
stations       <- sim$design$stations       # one row per station

sp_colours <- viridis(n_species, option = "viridis", end = 0.85)

# =============================================================================
# Row 1 — Surface density on 1 km X/Y grid
#
# Deterministic mean surface:
#   log(lambda) = mu_s + zbathy_pref(Z_bathy(X)) + 0  [no GP residual]
# Z_bathy is derived from the shelf-slope profile function using X.
# Sampled station locations overlaid as open circles.
# =============================================================================

bathy_mean_fn <- function(X_km) {
  abyssal <- 2500
  shelf   <- 80
  pmax(abyssal + (shelf - abyssal) / (1 + exp(-0.06 * (X_km - 180))), 10)
}

zbathy_pref_fn <- function(Z_bathy, mu_z, sd_z, amp) {
  raw <- dnorm(Z_bathy, mu_z, sd_z)
  raw <- raw / max(raw)
  amp * (raw - 0.5)
}

grid_1km <- expand.grid(
  X = seq(0, 300, by = 1),
  Y = seq(0, 400, by = 1)
) %>%
  mutate(Z_bathy = bathy_mean_fn(X))

for (s in seq_len(n_species)) {
  p    <- gp_params[[s]]
  z_mu <- zbathy_pref_fn(grid_1km$Z_bathy,
                         p$zbathy_pref_mu,
                         p$zbathy_pref_sd,
                         p$zbathy_pref_amp)
  grid_1km[[paste0("lambda_s", s)]] <- exp(p$mu + z_mu)
}

# Shared colour limits (log scale)
log_lam_all <- unlist(lapply(seq_len(n_species), function(s)
  log(grid_1km[[paste0("lambda_s", s)]] + 0.01)))
clim <- quantile(log_lam_all, c(0.01, 0.99))

# Surface samples (Z_sample == 0) for overlay
surf_samples <- samples %>%
  filter(Z_sample == 0) %>%
  mutate(
    lambda_s1 = lambda_true[sample_id, 1],
    lambda_s2 = lambda_true[sample_id, 2],
    lambda_s3 = lambda_true[sample_id, 3]
  )

p_row1 <- lapply(seq_len(n_species), function(s) {
  
  grid_df <- grid_1km %>%
    transmute(X, Y,
              log_lambda = log(.data[[paste0("lambda_s", s)]] + 0.01))
  
  surf_pts <- surf_samples %>%
    transmute(X, Y,
              log_lambda_obs = log(.data[[paste0("lambda_s", s)]] + 0.01))
  
  ggplot() +
    geom_raster(
      data = grid_df,
      aes(X, Y, fill = log_lambda)
    ) +
    geom_point(
      data   = surf_pts,
      aes(X, Y, colour = log_lambda_obs),
      shape  = 21,
      size   = 1.8,
      fill   = NA,
      stroke = 0.7
    ) +
    # Approximate shelf break (~200 m isobath) at X ≈ 180 km
    geom_vline(
      xintercept = 180,
      linetype   = "dashed",
      colour     = "white",
      linewidth  = 0.5,
      alpha      = 0.75
    ) +
    annotate(
      "text", x = 183, y = 388,
      label  = "shelf break",
      colour = "white",
      size   = 2.6,
      hjust  = 0,
      fontface = "italic"
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = clim,
      name   = "log(λ)",
      oob    = scales::squish
    ) +
    scale_colour_viridis_c(
      option = "viridis",
      limits = clim,
      guide  = "none",
      oob    = scales::squish
    ) +
    scale_x_continuous(breaks = seq(0, 300, 100), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 400, 100), expand = c(0, 0)) +
    coord_fixed(ratio = 1, xlim = c(0, 300), ylim = c(0, 400)) +
    labs(
      title    = sp_common[s],
      subtitle = "Mean surface density (Z_sample = 0 m)",
      x = "Easting (km from offshore)",
      y = "Northing (km from Cape Blanco)"
    ) +
    theme_bw(base_size = 11) +
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
# Row 2 — Bathymetry (Z_bathy) vs true animal density
# =============================================================================

# Build a long data frame directly from samples without the broken select()
bathy_df <- bind_rows(lapply(seq_len(n_species), function(s) {
  data.frame(
    station = samples$station,
    Z_bathy = samples$Z_bathy,
    lambda  = lambda_true[, s],
    species = sp_common[s]
  )
})) %>%
  # Average lambda over sample depths within each station
  group_by(station, Z_bathy, species) %>%
  summarise(lambda = mean(lambda), .groups = "drop") %>%
  mutate(species = factor(species, levels = sp_common))

# Deterministic mean function over a fine Z_bathy grid
zbathy_seq <- seq(10, 2600, by = 10)

bathy_mean_df <- bind_rows(lapply(seq_len(n_species), function(s) {
  p    <- gp_params[[s]]
  z_mu <- zbathy_pref_fn(zbathy_seq,
                         p$zbathy_pref_mu,
                         p$zbathy_pref_sd,
                         p$zbathy_pref_amp)
  data.frame(
    Z_bathy = zbathy_seq,
    lambda  = exp(p$mu + z_mu),
    species = sp_common[s]
  )
})) %>%
  mutate(species = factor(species, levels = sp_common))

p_row2 <- ggplot() +
  geom_point(
    data  = bathy_df,
    aes(x = Z_bathy, y = lambda, colour = species),
    size  = 1.2,
    alpha = 0.35
  ) +
  geom_line(
    data      = bathy_mean_df,
    aes(x = Z_bathy, y = lambda, colour = species),
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = 200,
    linetype   = "dashed",
    colour     = "grey50",
    linewidth  = 0.5
  ) +
  annotate(
    "text", x = 210, y = max(bathy_mean_df$lambda) * 0.95,
    label    = "shelf break\n(~200 m)",
    colour   = "grey40",
    size     = 3.0,
    hjust    = 0,
    fontface = "italic"
  ) +
  scale_colour_manual(
    values = sp_colours,
    name   = NULL,
    labels = sp_common
  ) +
  scale_x_continuous(
    breaks = c(0, 200, 500, 1000, 1500, 2000, 2500),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(expand = c(0.02, 0)) +
  labs(
    title    = "Bathymetric habitat preference",
    subtitle = paste("True animal density λ vs bottom depth Z_bathy;",
                     "line = deterministic mean function,",
                     "points = station means (including GP residual)"),
    x = "Bottom depth Z_bathy (m)",
    y = "λ (copies / L)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(size = 11, face = "bold"),
    plot.subtitle   = element_text(size = 9, colour = "grey40"),
    legend.position = "right",
    legend.text     = element_text(face = "italic", size = 9)
  )

# =============================================================================
# Row 3 — Water column sampling depth effect (Z_sample)
#
# Shows exp(zsample_pref[d]) for each species at each sample depth.
# Depth on reversed y-axis (oceanographic convention: shallow at top).
# This is entirely separate from animal density — it reflects where in
# the water column eDNA is concentrated given animals are present.
# =============================================================================

zsample_df <- bind_rows(lapply(seq_len(n_species), function(s) {
  p <- gp_params[[s]]
  tibble(
    species        = sp_common[s],
    depth_idx      = seq_len(D),
    Z_sample       = sample_depths,
    log_offset     = p$zsample_pref,
    edna_mult      = exp(p$zsample_pref)
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
    geom_vline(
      xintercept = 1,
      linetype   = "dashed",
      colour     = "grey60",
      linewidth  = 0.5
    ) +
    geom_path(
      colour      = col,
      linewidth   = 1.0,
      orientation = "y"
    ) +
    geom_point(colour = col, size = 3.5) +
    geom_text(
      aes(label = sprintf("%+.2f", log_offset)),
      hjust  = -0.30,
      size   = 3.2,
      colour = "grey25"
    ) +
    scale_y_reverse(
      breaks = sample_depths,
      labels = paste0(sample_depths, " m")
    ) +
    scale_x_continuous(
      limits = xlim_mult,
      expand = c(0.05, 0)
    ) +
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
  plot_layout(heights = c(1.4, 0.9, 1.0)) +
  plot_annotation(
    title    = "Simulated eDNA fields — Oregon / Washington coast",
    subtitle = paste0(
      "Row 1: mean surface density (1 km grid, Z_sample = 0 m)  |  ",
      "Row 2: true animal density λ vs bathymetry Z_bathy  |  ",
      "Row 3: eDNA water column concentration effect by sample depth"
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  )

ggsave("simulated_edna_fields.png", fig,
       width = 14, height = 15, dpi = 150)

cat("Saved simulated_edna_fields.png\n")
print(fig)
