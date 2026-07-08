# =============================================================================
# scripts/plotHarborPorpoiseeDNA.R
#
# Single-panel eDNA-only version of scripts/plotHarborPorpoiseData.R for
# harbor porpoise (Phocoena phocoena). Shows the MARVER1 metabarcoding
# sampled locations (dark crosses) and harbor porpoise detections (red
# circles sized by total target reads).
#
# Style mirrors the right panel of plotHarborPorpoiseData.R:
#   * Light-grey land (CA/OR/WA from `state` + Canada from `world`)
#   * Dark crosses  = sampled locations
#   * Red circles   = detections, sqrt-scaled point size by total reads
#
# Inputs : data/MV1_MURI_df.csv
# Output : figures/edna_harbor_porpoise.png  (8 x 8 in, 200 dpi)
#
# Run from the project root:
#   Rscript scripts/plotHarborPorpoiseeDNA.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(maps)
})

# ---------------------------------------------------------------------------
# Constants — species + map extent (San Francisco and northward)
# ---------------------------------------------------------------------------
SPECIES_LABEL <- "Harbor porpoise"
MB_SPECIES    <- "Phocoena phocoena"

LAT_MIN <- 37.5
LON_MIN <- -130
LON_MAX <- -120

OUT_PNG <- "figures/edna_harbor_porpoise.png"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
mv1 <- read.csv("data/MV1_MURI_df.csv", stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# eDNA subsets — all MARVER1 sampling locations + harbor porpoise detections
# aggregated to (location_id) with total reads.
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

cat(sprintf("eDNA: %d sampled locations, %d detections (max reads = %s)\n",
            nrow(mb_all_loc), nrow(mb_sp_loc),
            ifelse(nrow(mb_sp_loc) == 0, "0",
                   format(max(mb_sp_loc$total_reads), big.mark = ","))))

# ---------------------------------------------------------------------------
# Map upper-latitude bound from data (no artificial cap), with headroom
# ---------------------------------------------------------------------------
LAT_MAX <- max(c(mb_all_loc$lat, mb_sp_loc$lat), na.rm = TRUE) + 0.5

# ---------------------------------------------------------------------------
# Coastline — state polygons for CA/OR/WA + Canada from world.
# ---------------------------------------------------------------------------
states <- map_data("state") |>
  filter(region %in% c("california", "oregon", "washington"))
canada <- map_data("world") |>
  filter(region == "Canada")

pretty_size_breaks <- function(x, n = 4) {
  br <- pretty(x, n = n)
  br <- br[br > 0 & br <= max(x)]
  if (length(br) < 2) br <- range(x)
  br
}

mb_breaks <- pretty_size_breaks(mb_sp_loc$total_reads)

# ---------------------------------------------------------------------------
# Build the plot
# ---------------------------------------------------------------------------
fig <- ggplot() +
  geom_polygon(data = canada, aes(long, lat, group = group),
               fill = "grey92", colour = "grey55", linewidth = 0.3) +
  geom_polygon(data = states, aes(long, lat, group = group),
               fill = "grey92", colour = "grey55", linewidth = 0.3) +
  geom_point(data = mb_all_loc,
             aes(x = lon, y = lat),
             shape = 4, size = 1.0,
             colour = "grey25", alpha = 0.7) +
  geom_point(data = mb_sp_loc,
             aes(x = lon, y = lat, size = total_reads),
             colour = "#C7383A", fill = "#C7383A", alpha = 0.55,
             shape = 21, stroke = 0.4) +
  coord_quickmap(xlim = c(LON_MIN, LON_MAX),
                 ylim = c(LAT_MIN, LAT_MAX),
                 expand = FALSE) +
  scale_size_continuous(name   = "Total reads",
                        range  = c(1.0, 9),
                        trans  = "sqrt",
                        breaks = mb_breaks) +
  scale_x_continuous(breaks = seq(-130, -120, 2),
                     labels = function(x) sprintf("%g°W", abs(x))) +
  scale_y_continuous(breaks = seq(38, 50, 2),
                     labels = function(y) sprintf("%g°N", y)) +
  labs(title    = sprintf("%s — MARVER1 eDNA observations", SPECIES_LABEL),
       subtitle = sprintf("%d detections / %d sampled locations. Dark crosses = all sampled locations; red circles = detections (sqrt-scaled by total reads, 0 – %s).",
                          nrow(mb_sp_loc), nrow(mb_all_loc),
                          format(max(mb_sp_loc$total_reads), big.mark = ",")),
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9, colour = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position  = "right")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PNG, fig, width = 8, height = 8, dpi = 200)
cat(sprintf("Saved %s\n", OUT_PNG))
