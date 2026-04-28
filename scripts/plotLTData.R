# =============================================================================
# scripts/plotLTData.R
#
# Two-panel map of line-transect survey effort + cetacean sightings on the
# US West Coast (San Francisco to the US/Canada border, ~37.5 - 49.5 deg N).
# One panel per species, with sighting points sized by group_size and
# INDEPENDENT group-size scales per panel:
#   * Pacific white-sided dolphin (spcode = 22; group sizes 1 - 452 in
#     these data, so a sqrt-transformed point-size scale)
#   * Humpback whale               (spcode = 76; group sizes 1 - 8,
#     linear point-size scale)
#
# Inputs : Data/effort.csv, Data/sightings.csv
# Output : figures/lt_data_pwsd_humpback.png  (12 x 7 in, 200 dpi)
#
# Run from the project root:
#   Rscript scripts/plotLTData.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(maps)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Map extent: SF (~37.77 deg N) and northward; no fixed upper latitude limit
# (data extent drives the upper bound).
# ---------------------------------------------------------------------------
LAT_MIN <- 37.5
LON_MIN <- -130
LON_MAX <- -120

# ---------------------------------------------------------------------------
# Load + filter
# ---------------------------------------------------------------------------
sights <- read.csv("Data/sightings.csv", stringsAsFactors = FALSE)
eff    <- read.csv("Data/effort.csv",    stringsAsFactors = FALSE)

# Sightings: keep PWSD (22) and Humpback (76), inside the map extent, with
# a non-missing positive group_size.
sp_codes <- c("Pacific white-sided dolphin" = 22L,
              "Humpback whale"              = 76L)
sp_levels <- names(sp_codes)

sights_sub <- sights |>
  filter(spcode %in% sp_codes,
         !is.na(group_size), group_size > 0,
         !is.na(lat), !is.na(lon),
         lat >= LAT_MIN,
         lon >= LON_MIN, lon <= LON_MAX) |>
  mutate(species = factor(case_when(
    spcode == 22L ~ "Pacific white-sided dolphin",
    spcode == 76L ~ "Humpback whale"
  ), levels = sp_levels))

cat(sprintf("Filtered sightings: %d total\n", nrow(sights_sub)))
print(sights_sub |>
  group_by(species) |>
  summarise(n_sightings = dplyr::n(),
            min_size = min(group_size), median_size = median(group_size),
            mean_size = mean(group_size), max_size = max(group_size),
            .groups = "drop"))

# Effort: keep segments whose midpoint falls inside the map extent. Use
# midpoint columns where present; fall back to the begin point.
eff_sub <- eff |>
  mutate(mlat_use = coalesce(mlat, latitude_begin),
         mlon_use = coalesce(mlon, longitude_begin)) |>
  filter(!is.na(mlat_use), !is.na(mlon_use),
         mlat_use >= LAT_MIN,
         mlon_use >= LON_MIN, mlon_use <= LON_MAX)

# Compute upper-latitude bound from the data so the plot has no
# artificial cap. Add a half-degree of headroom for visual breathing
# room above the northernmost feature.
LAT_MAX <- max(c(sights_sub$lat, eff_sub$mlat_use), na.rm = TRUE) + 0.5

cat(sprintf("Filtered effort segments: %d (of %d total)\n",
            nrow(eff_sub), nrow(eff)))

# ---------------------------------------------------------------------------
# Coastline. Use `state` for the US (gives CA/OR/WA shoreline + state
# borders) and `world` filtered to Canada for the bit north of WA.
# Plotting both `world$USA` and `state$ca/or/wa` in earlier versions
# produced doubled outlines because both polygons traced the same coast.
# ---------------------------------------------------------------------------
states <- map_data("state") |>
  filter(region %in% c("california", "oregon", "washington"))
canada <- map_data("world") |>
  filter(region == "Canada")

# ---------------------------------------------------------------------------
# Per-species panel builder
# ---------------------------------------------------------------------------
make_panel <- function(species_label, sights_data, size_range, size_trans = "identity") {
  ggplot() +
    # Land: Canada (north of WA) + CA/OR/WA via state polygons.
    # Plotting only one source for each region avoids the doubled-
    # outline effect from earlier versions.
    geom_polygon(data = canada, aes(long, lat, group = group),
                 fill = "grey92", colour = "grey55", linewidth = 0.3) +
    geom_polygon(data = states, aes(long, lat, group = group),
                 fill = "grey92", colour = "grey55", linewidth = 0.3) +
    # Effort segments (begin->end). Darker / thicker than the first
    # version so they are visible against the land + ocean.
    geom_segment(data = eff_sub,
                 aes(x = longitude_begin, y = latitude_begin,
                     xend = longitude_end, yend = latitude_end),
                 colour = "grey25", linewidth = 0.55, alpha = 0.7) +
    # Sighting points sized by group_size
    geom_point(data = sights_data,
               aes(x = lon, y = lat, size = group_size),
               colour = "#C7383A", fill = "#C7383A", alpha = 0.55,
               shape = 21, stroke = 0.4) +
    coord_quickmap(xlim = c(LON_MIN, LON_MAX),
                   ylim = c(LAT_MIN, LAT_MAX),
                   expand = FALSE) +
    scale_size_continuous(name = "Group size",
                          range = size_range,
                          trans = size_trans,
                          breaks = NULL) +
    scale_x_continuous(breaks = seq(-130, -120, 2),
                       labels = function(x) sprintf("%g°W", abs(x))) +
    scale_y_continuous(breaks = seq(38, 49, 2),
                       labels = function(y) sprintf("%g°N", y)) +
    labs(title    = sprintf("%s (n = %d sightings)",
                            species_label, nrow(sights_data)),
         subtitle = sprintf("group size %d - %d",
                            min(sights_data$group_size),
                            max(sights_data$group_size)),
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 12),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank(),
          legend.position  = "right")
}

# Per-panel breaks (each panel will also auto-pick from data, so this
# provides nice round numbers for the legend)
make_panel_breaks <- function(species_label, sights_data, size_range, size_trans = "identity") {
  pretty_breaks <- pretty(sights_data$group_size, n = 4)
  pretty_breaks <- pretty_breaks[pretty_breaks > 0 & pretty_breaks <= max(sights_data$group_size)]
  if (length(pretty_breaks) < 2) pretty_breaks <- range(sights_data$group_size)

  make_panel(species_label, sights_data, size_range, size_trans) +
    scale_size_continuous(name = "Group size",
                          range = size_range,
                          trans = size_trans,
                          breaks = pretty_breaks)
}

# ---------------------------------------------------------------------------
# Build panels with INDEPENDENT group-size scales:
#   PWSD has very wide group-size range (1-452). Use sqrt-transformed sizes
#     so the small groups are still visible.
#   Humpback group sizes are 1-8 (mostly 1s). Linear scale is fine.
# ---------------------------------------------------------------------------
sights_pwsd     <- filter(sights_sub, species == "Pacific white-sided dolphin")
sights_humpback <- filter(sights_sub, species == "Humpback whale")

p_pwsd <- make_panel_breaks("Pacific white-sided dolphin",
                            sights_pwsd,
                            size_range = c(0.6, 9),
                            size_trans = "sqrt")
p_humpback <- make_panel_breaks("Humpback whale",
                                sights_humpback,
                                size_range = c(1.2, 6),
                                size_trans = "identity")

# ---------------------------------------------------------------------------
# Combine
# ---------------------------------------------------------------------------
fig <- p_pwsd + p_humpback +
  plot_annotation(
    title    = "Line-transect effort + cetacean sightings",
    subtitle = "US West Coast: San Francisco to US/Canada border. Grey lines = effort segments; red points = sightings (sized by group size, independent scales per panel).",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10, colour = "grey40"))
  )

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
out_png <- "figures/lt_data_pwsd_humpback.png"
ggsave(out_png, fig, width = 12, height = 7, dpi = 200)
cat(sprintf("Saved %s\n", out_png))
