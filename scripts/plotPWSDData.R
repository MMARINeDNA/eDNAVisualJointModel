# =============================================================================
# scripts/plotPWSDData.R
#
# Two-panel side-by-side figure for Pacific white-sided dolphin
# (Lagenorhynchus obliquidens), combining:
#
#   Left panel  — Line-transect effort + PWSD sightings (sized by
#                 group_size; sqrt-transformed because group sizes
#                 span 1 – 452).
#   Right panel — MARVER1 metabarcoding sampling locations + PWSD
#                 detections (sized by total target reads;
#                 sqrt-transformed because reads span 0 – ~50k).
#
# Style mirrors scripts/plotLTData.R + scripts/ploteDNAData.R:
#   * Light-grey land (CA/OR/WA from `state` + Canada from `world`)
#   * LT panel:    grey segments = effort, red sized circles = sightings
#   * eDNA panel:  dark crosses  = sampled locations, red sized circles
#                                  = detections
#   * Independent legends per panel (units differ: group size vs reads)
#
# Inputs : data/effort.csv, data/sightings.csv, data/MV1_MURI_df.csv
# Output : figures/lt_edna_pwsd.png  (14 x 7 in, 200 dpi)
#
# Run from the project root:
#   Rscript scripts/plotPWSDData.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(maps)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Constants — species + map extent (San Francisco and northward)
# ---------------------------------------------------------------------------
SPECIES_LABEL  <- "Pacific white-sided dolphin"
LT_SPCODE      <- 22L
MB_SPECIES     <- "Lagenorhynchus obliquidens"

LAT_MIN <- 37.5
LON_MIN <- -130
LON_MAX <- -120

OUT_PNG <- "figures/lt_edna_pwsd.png"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
sights <- read.csv("data/sightings.csv",      stringsAsFactors = FALSE)
eff    <- read.csv("data/effort.csv",         stringsAsFactors = FALSE)
mv1    <- read.csv("data/MV1_MURI_df.csv",    stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# LT subsets — PWSD sightings + all effort segments inside the map extent
# ---------------------------------------------------------------------------
sights_sp <- sights |>
  filter(spcode == LT_SPCODE,
         !is.na(group_size), group_size > 0,
         !is.na(lat), !is.na(lon),
         lat >= LAT_MIN,
         lon >= LON_MIN, lon <= LON_MAX)

eff_sub <- eff |>
  mutate(mlat_use = coalesce(mlat, latitude_begin),
         mlon_use = coalesce(mlon, longitude_begin)) |>
  filter(!is.na(mlat_use), !is.na(mlon_use),
         mlat_use >= LAT_MIN,
         mlon_use >= LON_MIN, mlon_use <= LON_MAX)

cat(sprintf("LT  : %d sightings, %d effort segments\n",
            nrow(sights_sp), nrow(eff_sub)))

# ---------------------------------------------------------------------------
# eDNA subsets — all MARVER1 sampling locations + PWSD detections aggregated
# to (location_id) with total reads.
# ---------------------------------------------------------------------------
mb_all_loc <- mv1 |>
  filter(!is.na(lat), !is.na(lon),
         lat >= LAT_MIN,
         lon >= LON_MIN, lon <= LON_MAX) |>
  distinct(location_id, lat, lon)

mb_sp_loc <- mv1 |>
  filter(species == MB_SPECIES,
         !is.na(lat), !is.na(lon),
         lat >= LAT_MIN,
         lon >= LON_MIN, lon <= LON_MAX) |>
  group_by(location_id) |>
  summarise(lat         = first(lat),
            lon         = first(lon),
            total_reads = sum(Nreads, na.rm = TRUE),
            .groups     = "drop") |>
  filter(total_reads > 0)

cat(sprintf("eDNA: %d sampled locations, %d detections (max reads = %d)\n",
            nrow(mb_all_loc), nrow(mb_sp_loc),
            ifelse(nrow(mb_sp_loc) == 0, 0L, max(mb_sp_loc$total_reads))))

# ---------------------------------------------------------------------------
# Map upper-latitude bound from data (no artificial cap), with headroom
# ---------------------------------------------------------------------------
LAT_MAX <- max(c(sights_sp$lat, eff_sub$mlat_use,
                 mb_all_loc$lat, mb_sp_loc$lat),
               na.rm = TRUE) + 0.5

# ---------------------------------------------------------------------------
# Coastline — state polygons for CA/OR/WA + Canada from world. Plotting
# only one source per region avoids the doubled-outline effect.
# ---------------------------------------------------------------------------
states <- map_data("state") |>
  filter(region %in% c("california", "oregon", "washington"))
canada <- map_data("world") |>
  filter(region == "Canada")

# Common ggplot scaffolding (land + coords + axes + theme) shared by both
# panels. Each panel adds its own background-effort/sample layer, its own
# detections layer, and its own size legend.
base_layers <- list(
  geom_polygon(data = canada, aes(long, lat, group = group),
               fill = "grey92", colour = "grey55", linewidth = 0.3),
  geom_polygon(data = states, aes(long, lat, group = group),
               fill = "grey92", colour = "grey55", linewidth = 0.3),
  coord_quickmap(xlim = c(LON_MIN, LON_MAX),
                 ylim = c(LAT_MIN, LAT_MAX),
                 expand = FALSE),
  scale_x_continuous(breaks = seq(-130, -120, 2),
                     labels = function(x) sprintf("%g°W", abs(x))),
  scale_y_continuous(breaks = seq(38, 50, 2),
                     labels = function(y) sprintf("%g°N", y)),
  labs(x = NULL, y = NULL),
  theme_bw(base_size = 11),
  theme(plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 9, colour = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position  = "right")
)

pretty_size_breaks <- function(x, n = 4) {
  br <- pretty(x, n = n)
  br <- br[br > 0 & br <= max(x)]
  if (length(br) < 2) br <- range(x)
  br
}

# ---------------------------------------------------------------------------
# Left panel — LT effort + sightings
# ---------------------------------------------------------------------------
lt_breaks <- pretty_size_breaks(sights_sp$group_size)

p_lt <- ggplot() +
  base_layers +
  geom_segment(data = eff_sub,
               aes(x = longitude_begin, y = latitude_begin,
                   xend = longitude_end, yend = latitude_end),
               colour = "grey25", linewidth = 0.55, alpha = 0.7) +
  geom_point(data = sights_sp,
             aes(x = lon, y = lat, size = group_size),
             colour = "#C7383A", fill = "#C7383A", alpha = 0.55,
             shape = 21, stroke = 0.4) +
  scale_size_continuous(name   = "Group size",
                        range  = c(0.6, 9),
                        trans  = "sqrt",
                        breaks = lt_breaks) +
  labs(title    = sprintf("Line-transect: %d sightings", nrow(sights_sp)),
       subtitle = sprintf("group size %d – %d (sqrt-scaled point size)",
                          min(sights_sp$group_size),
                          max(sights_sp$group_size)))

# ---------------------------------------------------------------------------
# Right panel — eDNA sampling + detections
# ---------------------------------------------------------------------------
mb_breaks <- pretty_size_breaks(mb_sp_loc$total_reads)

p_edna <- ggplot() +
  base_layers +
  geom_point(data = mb_all_loc,
             aes(x = lon, y = lat),
             shape = 4, size = 1.0,
             colour = "grey25", alpha = 0.7) +
  geom_point(data = mb_sp_loc,
             aes(x = lon, y = lat, size = total_reads),
             colour = "#C7383A", fill = "#C7383A", alpha = 0.55,
             shape = 21, stroke = 0.4) +
  scale_size_continuous(name   = "Total reads",
                        range  = c(1.0, 9),
                        trans  = "sqrt",
                        breaks = mb_breaks) +
  labs(title    = sprintf("MARVER1 eDNA: %d detections / %d sampled locations",
                          nrow(mb_sp_loc), nrow(mb_all_loc)),
       subtitle = sprintf("total reads 0 – %s (sqrt-scaled point size)",
                          format(max(mb_sp_loc$total_reads), big.mark = ",")))

# ---------------------------------------------------------------------------
# Combine
# ---------------------------------------------------------------------------
fig <- p_lt + p_edna +
  plot_annotation(
    title    = sprintf("%s line-transect + MARVER1 eDNA observations",
                       SPECIES_LABEL),
    subtitle = "Left: grey lines = LT effort, red points = sightings. Right: dark crosses = eDNA sampled locations, red points = detections.",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10, colour = "grey40"))
  )

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PNG, fig, width = 14, height = 7, dpi = 200)
cat(sprintf("Saved %s\n", OUT_PNG))
