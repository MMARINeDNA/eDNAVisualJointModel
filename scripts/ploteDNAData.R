# =============================================================================
# scripts/ploteDNAData.R
#
# Four-panel map of eDNA observations on the US West Coast (San Francisco
# and northward), in the same style as scripts/plotLTData.R:
#
#   1. Pacific hake qPCR  (Data/hake_qPCR_MURI_df.csv)
#   2. Pacific hake               via metabarcoding (Data/MV1_MURI_df.csv)
#   3. Pacific white-sided dolphin via metabarcoding
#   4. Humpback whale              via metabarcoding
#
# Each panel:
#   * Light-grey land (CA/OR/WA states + Canada)
#   * Sampling locations shown as small dark crosses (analog of effort
#     segments in the LT plot - "where we sampled")
#   * Detections shown as red filled circles, SIZED by an abundance
#     measure with a panel-specific scale (panel-independent legends):
#       qPCR hake: number of positive qPCR replicates at the location
#       MB species: total target reads at the location (sqrt-transformed
#                   because the range spans 0 to ~10^6)
#
# Inputs : Data/hake_qPCR_MURI_df.csv, Data/MV1_MURI_df.csv
# Output : figures/edna_data_4panel.png  (14 x 10 in, 200 dpi)
#
# Run from the project root:
#   Rscript scripts/ploteDNAData.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(maps)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Map extent (matches plotLTData.R: SF and northward, no upper cap)
# ---------------------------------------------------------------------------
LAT_MIN <- 37.5
LON_MIN <- -130
LON_MAX <- -120

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
qpcr <- read.csv("Data/hake_qPCR_MURI_df.csv", stringsAsFactors = FALSE)
mv1  <- read.csv("Data/MV1_MURI_df.csv",       stringsAsFactors = FALSE)

cat(sprintf("hake_qPCR rows: %d\n", nrow(qpcr)))
cat(sprintf("MV1 rows      : %d\n", nrow(mv1)))

# ---------------------------------------------------------------------------
# Aggregate qPCR hake: one row per location_id with detection counts
# ---------------------------------------------------------------------------
qpcr_loc <- qpcr |>
  filter(!is.na(lat), !is.na(lon),
         lat >= LAT_MIN,
         lon >= LON_MIN, lon <= LON_MAX) |>
  group_by(location_id) |>
  summarise(lat       = first(lat),
            lon       = first(lon),
            n_reps    = dplyr::n(),
            n_pos     = sum(detected, na.rm = TRUE),
            max_conc  = max(conc,     na.rm = TRUE),
            .groups   = "drop") |>
  mutate(detected_any = n_pos > 0)

cat(sprintf("qPCR locations: %d (%d with at least one positive)\n",
            nrow(qpcr_loc), sum(qpcr_loc$detected_any)))

# ---------------------------------------------------------------------------
# Aggregate MV1: one row per (location_id, species) with total target reads
# ---------------------------------------------------------------------------
mb_targets <- c("Merluccius productus",       # hake
                "Lagenorhynchus obliquidens", # PWSD
                "Megaptera novaeangliae")     # humpback

mb_loc <- mv1 |>
  filter(species %in% mb_targets,
         !is.na(lat), !is.na(lon),
         lat >= LAT_MIN,
         lon >= LON_MIN, lon <= LON_MAX) |>
  group_by(species, location_id) |>
  summarise(lat        = first(lat),
            lon        = first(lon),
            total_reads = sum(Nreads, na.rm = TRUE),
            .groups   = "drop") |>
  mutate(detected_any = total_reads > 0)

# Common set of MB sampling locations (across all target species these
# are the same 551 locations in MV1; collapse for the "all sampled
# locations" overlay).
mb_all_loc <- mv1 |>
  filter(!is.na(lat), !is.na(lon),
         lat >= LAT_MIN,
         lon >= LON_MIN, lon <= LON_MAX) |>
  distinct(location_id, lat, lon)

cat(sprintf("MB locations: %d total (across all 295 species)\n",
            nrow(mb_all_loc)))
mb_loc |>
  group_by(species) |>
  summarise(n_locations = dplyr::n(),
            n_detected  = sum(detected_any),
            max_reads   = max(total_reads),
            .groups     = "drop") |>
  print()

# ---------------------------------------------------------------------------
# Compute LAT_MAX from data so the plot has no artificial upper cap
# ---------------------------------------------------------------------------
LAT_MAX <- max(c(qpcr_loc$lat, mb_all_loc$lat), na.rm = TRUE) + 0.5

# ---------------------------------------------------------------------------
# Coastline. Use `state` for the US (CA/OR/WA) and `world` for Canada;
# this avoids the doubled-outline effect that came from plotting both
# `world$USA` and the state polygons.
# ---------------------------------------------------------------------------
states <- map_data("state") |>
  filter(region %in% c("california", "oregon", "washington"))
canada <- map_data("world") |>
  filter(region == "Canada")

# ---------------------------------------------------------------------------
# Per-panel builder. `bg_data` is the "all sampling locations" overlay
# (small dark crosses); `det_data` is the per-detection points sized by
# `size_var`.
# ---------------------------------------------------------------------------
make_panel <- function(panel_label, bg_data, det_data, size_var,
                       size_label, size_range, size_trans = "identity") {
  pretty_breaks <- pretty(det_data[[size_var]], n = 4)
  pretty_breaks <- pretty_breaks[pretty_breaks > 0 &
                                 pretty_breaks <= max(det_data[[size_var]])]
  if (length(pretty_breaks) < 2) {
    pretty_breaks <- range(det_data[[size_var]])
  }

  ggplot() +
    # Land
    geom_polygon(data = canada, aes(long, lat, group = group),
                 fill = "grey92", colour = "grey55", linewidth = 0.3) +
    geom_polygon(data = states, aes(long, lat, group = group),
                 fill = "grey92", colour = "grey55", linewidth = 0.3) +
    # All sampling locations: small dark crosses (analog of effort)
    geom_point(data = bg_data,
               aes(x = lon, y = lat),
               shape = 4, size = 1.0,
               colour = "grey25", alpha = 0.7) +
    # Detections: sized red points
    geom_point(data = det_data,
               aes(x = lon, y = lat, size = .data[[size_var]]),
               colour = "#C7383A", fill = "#C7383A", alpha = 0.55,
               shape = 21, stroke = 0.4) +
    coord_quickmap(xlim = c(LON_MIN, LON_MAX),
                   ylim = c(LAT_MIN, LAT_MAX),
                   expand = FALSE) +
    scale_size_continuous(name   = size_label,
                          range  = size_range,
                          trans  = size_trans,
                          breaks = pretty_breaks) +
    scale_x_continuous(breaks = seq(-130, -120, 2),
                       labels = function(x) sprintf("%g°W", abs(x))) +
    scale_y_continuous(breaks = seq(38, 50, 2),
                       labels = function(y) sprintf("%g°N", y)) +
    labs(title    = panel_label,
         subtitle = sprintf("%d sampled locations, %d positive (detections shown as red points)",
                            nrow(bg_data), nrow(det_data)),
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank(),
          legend.position  = "right")
}

# ---------------------------------------------------------------------------
# Panel data subsets
# ---------------------------------------------------------------------------
qpcr_det <- filter(qpcr_loc, detected_any)
qpcr_bg  <- qpcr_loc                              # all qPCR locations

mb_hake_det     <- filter(mb_loc,
                          species == "Merluccius productus", detected_any)
mb_pwsd_det     <- filter(mb_loc,
                          species == "Lagenorhynchus obliquidens", detected_any)
mb_humpback_det <- filter(mb_loc,
                          species == "Megaptera novaeangliae", detected_any)
mb_bg           <- mb_all_loc                     # all MB locations

# ---------------------------------------------------------------------------
# Build panels (independent legends per panel; each gets its own
# size_range / transform tuned to the abundance range of that species)
# ---------------------------------------------------------------------------
p1 <- make_panel(
  panel_label = "Pacific hake — qPCR",
  bg_data    = qpcr_bg,
  det_data   = qpcr_det,
  size_var   = "n_pos",
  size_label = "qPCR positives",
  size_range = c(1.5, 7.5),
  size_trans = "identity"
)

p2 <- make_panel(
  panel_label = "Pacific hake — MARVER",
  bg_data    = mb_bg,
  det_data   = mb_hake_det,
  size_var   = "total_reads",
  size_label = "Total reads",
  size_range = c(0.8, 9),
  size_trans = "sqrt"
)

p3 <- make_panel(
  panel_label = "Pacific white-sided dolphin — MARVER",
  bg_data    = mb_bg,
  det_data   = mb_pwsd_det,
  size_var   = "total_reads",
  size_label = "Total reads",
  size_range = c(0.8, 9),
  size_trans = "sqrt"
)

p4 <- make_panel(
  panel_label = "Humpback whale — MARVER",
  bg_data    = mb_bg,
  det_data   = mb_humpback_det,
  size_var   = "total_reads",
  size_label = "Total reads",
  size_range = c(0.8, 9),
  size_trans = "sqrt"
)

# ---------------------------------------------------------------------------
# Combine into a 2 x 2 grid
# ---------------------------------------------------------------------------
fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title    = "eDNA observations: hake qPCR + three-species MARVER metabarcoding",
    subtitle = "US West Coast: San Francisco and northward. Dark crosses = all sampled locations; red circles = detections, sized by abundance (independent panel scales).",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10, colour = "grey40"))
  )

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
out_png <- "figures/edna_data_4panel.png"
ggsave(out_png, fig, width = 16, height = 11, dpi = 200)
cat(sprintf("Saved %s\n", out_png))
