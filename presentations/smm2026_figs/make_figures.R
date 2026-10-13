# =============================================================================
# presentations/smm2026_figs/make_figures.R
#
# Figures for presentations/smm2026.qmd (SMM 2026 talk).
#
#   1. Real-data maps on the 35°N study domain, from the committed processed
#      objects (data/processed/muri_{edna,visual}.rds), so the counts match
#      what the models are fitted to:
#        data_footprint.png   LT effort (2018) and eDNA bottles (2019)
#        data_humpback.png    humpback: LT sightings | metabarcoding detections
#        data_pwsd.png        PWSD:     LT sightings | metabarcoding detections
#        data_hake.png        hake:     qPCR detections | metabarcoding reads
#   2. hsgp_illustration.png: what an HSGP is (1-D basis sum; 2-D anisotropic
#      draws on the study box with the simulation length-scales).
#   3. PNG copies of the model-structure diagrams in docs/figures/, with the
#      prefers-color-scheme dark block removed (it would put dark boxes on a
#      white slide). Rasterised at 3x with headless Chrome for its font
#      fallback; rsvg is used if Chrome is missing, but drops the Unicode
#      subscripts.
#
# Run from the repo root:
#   Rscript presentations/smm2026_figs/make_figures.R
# =============================================================================

# Labels use non-ASCII characters (degree sign, en dash, script l); without a
# UTF-8 locale R replaces them with "..." on output.
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(maps)
  library(patchwork)
})

OUT <- "presentations/smm2026_figs"

COL_VIS  <- "#D55E00"   # visual sightings (Okabe-Ito vermillion)
COL_EDNA <- "#0072B2"   # eDNA detections  (Okabe-Ito blue)
COL_EFF  <- "grey55"
LAT_MIN  <- 35
XLIM <- c(-131.5, -121.5)
YLIM <- c(34.6, 50.8)

edna   <- readRDS("data/processed/muri_edna.rds")
visual <- readRDS("data/processed/muri_visual.rds")

# ---- Land -------------------------------------------------------------------
land <- bind_rows(
  map_data("world", region = c("USA", "Canada", "Mexico")),
  map_data("state", region = c("california", "oregon", "washington")) |>
    mutate(group = group + 1e5)
)

base_map <- function() {
  list(
    geom_polygon(data = land, aes(long, lat, group = group),
                 fill = "grey92", colour = "grey60", linewidth = 0.25),
    geom_hline(yintercept = LAT_MIN, linetype = "dashed", colour = "grey40",
               linewidth = 0.3),
    coord_quickmap(xlim = XLIM, ylim = YLIM, expand = FALSE),
    scale_x_continuous(breaks = seq(-130, -122, 4),
                       labels = function(x) paste0(abs(x), "\u00b0W")),
    scale_y_continuous(breaks = seq(36, 50, 2),
                       labels = function(x) paste0(x, "\u00b0N")),
    labs(x = NULL, y = NULL),
    theme_bw(base_size = 12),
    theme(panel.grid = element_line(colour = "grey94", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey30", size = 10.5),
          legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  )
}

# ---- Data tables ------------------------------------------------------------
segs <- visual$design$segments
bottles <- edna$design$samples

# Metabarcoding: total target reads per bottle (summed over PCR replicates)
mb <- as_tibble(edna$observed$mb_reads) |>
  mutate(sample_idx = edna$observed$mb_sample_idx) |>
  group_by(sample_idx) |>
  summarise(across(c(hake, humpback, pwsd), sum), .groups = "drop")
mb_bottles <- bottles |>
  mutate(sample_idx = row_number()) |>
  left_join(mb, by = "sample_idx")

# qPCR (hake): share of replicates positive per bottle
qpcr <- tibble(sample_idx = edna$observed$qpcr_sample_idx,
               det = edna$observed$qpcr_detect) |>
  group_by(sample_idx) |>
  summarise(n_rep = n(), n_pos = sum(det), .groups = "drop")
qpcr_bottles <- bottles |>
  mutate(sample_idx = row_number()) |>
  left_join(qpcr, by = "sample_idx")

# Location-level collapse for the eDNA maps: bottles at several depths share a
# station position, so plot one symbol per position with the depths pooled.
edna_locs <- function(df, value) {
  df |>
    group_by(lat = round(lat, 3), lon = round(lon, 3)) |>
    summarise(value = sum({{ value }}, na.rm = TRUE), .groups = "drop")
}

sightings <- function(sp) {
  visual$observed[[sp]]$obs |>
    filter(distance <= visual$meta$w[[sp]])
}

# ---- 1. Footprint -----------------------------------------------------------
n_locs <- nrow(distinct(bottles, round(lat, 3), round(lon, 3)))

p_foot <- ggplot() +
  base_map() +
  geom_point(data = segs, aes(lon, lat, colour = "Line-transect effort (2018)"),
             size = 0.35, shape = 16) +
  geom_point(data = distinct(bottles, lat, lon),
             aes(lon, lat, colour = "eDNA stations (2019)"),
             size = 1.2, shape = 4, stroke = 0.6) +
  scale_colour_manual(NULL,
                      values = c("Line-transect effort (2018)" = COL_EFF,
                                 "eDNA stations (2019)" = COL_EDNA),
                      breaks = c("Line-transect effort (2018)",
                                 "eDNA stations (2019)")) +
  guides(colour = guide_legend(override.aes = list(size = c(1.6, 2),
                                                   shape = c(16, 4)))) +
  labs(title = "Two surveys, one coast",
       subtitle = sprintf(
         "%s km of visual effort in %s segments\n%d eDNA bottles at %d stations, 0\u2013500 m depth",
         format(round(sum(segs$seg_l)), big.mark = ","),
         format(nrow(segs), big.mark = ","), nrow(bottles), n_locs))

ggsave(file.path(OUT, "data_footprint.png"), p_foot,
       width = 5.2, height = 7.8, dpi = 220, device = ragg::agg_png)

# ---- 2. Per-species pairs ---------------------------------------------------
species_pair <- function(sp, label) {
  obs <- sightings(sp)
  locs <- edna_locs(mb_bottles, .data[[sp]])
  n_bottle_pos <- sum(mb_bottles[[sp]] > 0, na.rm = TRUE)

  p_lt <- ggplot() +
    base_map() +
    geom_point(data = segs, aes(lon, lat), colour = COL_EFF, size = 0.25,
               shape = 16) +
    geom_point(data = obs, aes(lon, lat, size = size), colour = COL_VIS,
               alpha = 0.55, shape = 16) +
    scale_size_area("Group size", max_size = if (sp == "pwsd") 7 else 4.5,
                    breaks = if (sp == "pwsd") c(10, 50, 150) else c(1, 3, 6)) +
    labs(title = "Visual line transects",
         subtitle = sprintf("%d groups detected within %.1f km of the line",
                            nrow(obs), visual$meta$w[[sp]]))

  p_ed <- ggplot() +
    base_map() +
    geom_point(data = locs, aes(lon, lat), colour = "grey35", size = 0.9,
               shape = 4, stroke = 0.4) +
    geom_point(data = filter(locs, value > 0), aes(lon, lat, size = value),
               colour = COL_EDNA, alpha = 0.55, shape = 16) +
    scale_size_area("Reads", max_size = 7, labels = scales::label_comma(),
                    breaks = scales::breaks_pretty(3)) +
    labs(title = "eDNA metabarcoding (MarVer1)",
         subtitle = sprintf("Reads in %d of %d bottles",
                            n_bottle_pos, nrow(mb_bottles)))

  (p_lt | p_ed) + plot_annotation(
    title = label,
    theme = theme(plot.title = element_text(face = "bold.italic", size = 15)))
}

ggsave(file.path(OUT, "data_humpback.png"),
       species_pair("humpback", "Humpback whale"),
       width = 9.6, height = 7.8, dpi = 220, device = ragg::agg_png)
ggsave(file.path(OUT, "data_pwsd.png"),
       species_pair("pwsd", "Pacific white-sided dolphin"),
       width = 9.6, height = 7.8, dpi = 220, device = ragg::agg_png)

# ---- 3. Hake ----------------------------------------------------------------
q_locs <- qpcr_bottles |>
  group_by(lat = round(lat, 3), lon = round(lon, 3)) |>
  summarise(frac = sum(n_pos, na.rm = TRUE) / sum(n_rep, na.rm = TRUE),
            .groups = "drop")
hake_mb <- edna_locs(mb_bottles, hake)

p_hq <- ggplot() +
  base_map() +
  geom_point(data = q_locs, aes(lon, lat, colour = frac), size = 1.6,
             shape = 16) +
  scale_colour_viridis_c("Share of replicates positive", limits = c(0, 1),
                         option = "mako", direction = -1, end = 0.9,
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = 8, barheight = 0.6)) +
  labs(title = "qPCR",
       subtitle = sprintf("Detected in %d of %d bottles",
                          sum(qpcr_bottles$n_pos > 0, na.rm = TRUE),
                          sum(!is.na(qpcr_bottles$n_rep))))
p_hm <- ggplot() +
  base_map() +
  geom_point(data = hake_mb, aes(lon, lat), colour = "grey35", size = 0.9,
             shape = 4, stroke = 0.4) +
  geom_point(data = filter(hake_mb, value > 0), aes(lon, lat, size = value),
             colour = COL_EDNA, alpha = 0.45, shape = 16) +
  scale_size_area("Reads", max_size = 6, labels = scales::label_comma(),
                  breaks = scales::breaks_pretty(3)) +
  labs(title = "Metabarcoding (MarVer1)",
       subtitle = sprintf("Reads in %d of %d bottles",
                          sum(mb_bottles$hake > 0, na.rm = TRUE),
                          nrow(mb_bottles)))
ggsave(file.path(OUT, "data_hake.png"),
       (p_hq | p_hm) + plot_annotation(
         title = "Pacific hake",
         theme = theme(plot.title = element_text(face = "bold.italic", size = 15))),
       width = 9.6, height = 7.8, dpi = 220, device = ragg::agg_png)

# ---- 4. HSGP illustration ---------------------------------------------------
# Left: a 1-D HSGP draw as a weighted sum of sine basis functions. Right: draws
# of the 2-D anisotropic field on the study box with the simulation
# length-scales (50 km E-W, 150 km N-S). Standalone: same basis and spectral
# density as the models, but not the HSGP4eDNA code.
se_spd <- function(w, l, sigma = 1) sigma^2 * sqrt(2 * pi) * l * exp(-0.5 * (l * w)^2)
hsgp_phi <- function(x, L, j) sqrt(1 / L) * sin(pi * j * (x + L) / (2 * L))

set.seed(7)
y_km <- seq(0, 1770, length.out = 600)
Ly <- 1.5 * 1770 / 2
yc <- y_km - 1770 / 2
J1 <- 1:40
beta1 <- rnorm(length(J1))
terms <- sapply(J1, function(j)
  sqrt(se_spd(pi * j / (2 * Ly), 150)) * beta1[j] * hsgp_phi(yc, Ly, j))
f1 <- rowSums(terms)
show <- order(-apply(abs(terms), 2, max))[1:6]
d1 <- bind_rows(
  as_tibble(terms[, show], .name_repair = ~ paste0("b", show)) |>
    mutate(y = y_km) |>
    pivot_longer(-y, names_to = "j", values_to = "v") |>
    mutate(what = "basis"),
  tibble(y = y_km, v = f1, j = "sum", what = "sum")
)
p_1d <- ggplot(d1, aes(y, v, group = j)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_line(data = filter(d1, what == "basis"), colour = "grey60",
            linewidth = 0.5) +
  geom_line(data = filter(d1, what == "sum"), colour = "#534AB7",
            linewidth = 1.3) +
  labs(x = "Along-coast distance (km)", y = "f(y)",
       title = "A smooth surface as a sum of fixed waves",
       subtitle = "Grey: six of the weighted sine basis functions.\nPurple: the sum of all 40") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey30", size = 10.5))

Lx2 <- 1.5 * 690 / 2; Ly2 <- 1.5 * 1770 / 2
gx <- seq(0, 690, by = 10); gy <- seq(0, 1770, by = 10)
Mx <- 1:36; My <- 1:44
phx <- sapply(Mx, function(j) hsgp_phi(gx - 690 / 2, Lx2, j))
phy <- sapply(My, function(j) hsgp_phi(gy - 1770 / 2, Ly2, j))
spd2 <- outer(se_spd(pi * Mx / (2 * Lx2), 50), se_spd(pi * My / (2 * Ly2), 150))
d2 <- map_dfr(1:3, function(k) {
  B <- matrix(rnorm(length(Mx) * length(My)), length(Mx)) * sqrt(spd2)
  f <- phx %*% B %*% t(phy)                        # grid x by grid y
  expand_grid(Y = gy, X = gx) |>
    mutate(f = as.vector(f), draw = paste("Draw", k))   # x varies fastest
})
p_2d <- ggplot(d2, aes(X, Y, fill = f)) +
  geom_raster() +
  facet_wrap(~draw, nrow = 1) +
  scale_fill_viridis_c("f(X, Y)", option = "viridis") +
  coord_fixed(expand = FALSE) +
  scale_x_continuous(breaks = c(0, 300, 600)) +
  labs(x = "East–west (km)", y = "North–south (km)",
       title = "2-D field on the study box",
       subtitle = "ℓ = 50 km east–west, 150 km north–south") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey30", size = 10.5),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")
ggsave(file.path(OUT, "hsgp_illustration.png"),
       (p_1d | p_2d) + plot_layout(widths = c(1.1, 1)),
       width = 11, height = 6.2, dpi = 220, device = ragg::agg_png)

# ---- 5. Model diagrams ------------------------------------------------------
chrome <- "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
for (f in c("edna", "visual", "joint")) {
  svg <- readLines(sprintf("docs/figures/model_%s.svg", f), warn = FALSE)
  txt <- paste(svg, collapse = "\n")
  txt <- sub("(?s)\\s*@media \\(prefers-color-scheme: dark\\) \\{.*?\\n  \\}\\n",
             "\n", txt, perl = TRUE)
  out_svg <- file.path(tempdir(), sprintf("model_%s.svg", f))
  writeLines(txt, out_svg)

  wh <- as.integer(regmatches(txt, regexec('width="(\\d+)" height="(\\d+)"', txt))[[1]][2:3])
  out_png <- file.path(OUT, sprintf("model_%s.png", f))
  if (file.exists(chrome)) {
    system2(chrome, c("--headless", "--disable-gpu", "--hide-scrollbars",
                      "--force-device-scale-factor=3",
                      sprintf("--window-size=%d,%d", wh[1], wh[2]),
                      sprintf("--screenshot=%s", normalizePath(out_png, mustWork = FALSE)),
                      paste0("file://", normalizePath(out_svg))),
            stdout = FALSE, stderr = FALSE)
  } else {
    rsvg::rsvg_png(out_svg, out_png, width = 3 * wh[1])
  }
}

# ---- Trim white margins (patchwork + fixed-aspect maps leave wide ones) -------
if (nzchar(Sys.which("magick"))) {
  for (png in list.files(OUT, "\\.png$", full.names = TRUE))
    system2("magick", c(shQuote(png), "-trim", "+repage", "-bordercolor", "white",
                        "-border", "30", shQuote(png)))
}

cat("Figures written to", OUT, "\n")
