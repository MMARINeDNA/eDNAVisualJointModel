# =============================================================================
# presentations/smm2026_figs/make_results.R
#
# Simulation-results figure for presentations/smm2026.qmd, from one joint
# comparison (R/compare_joint.r; outputs are git-ignored, so run this where
# they exist):
#   results_field.png   posterior mean against true log density at every
#                       eDNA station and transect segment location, species
#                       (rows) x fit (columns), with field R2 and mean 95%
#                       interval width per panel
#
# Run from the repo root:
#   Rscript presentations/smm2026_figs/make_results.R
# Env: CMP_DIR (default outputs/joint/compare_s1)
# =============================================================================

invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
suppressPackageStartupMessages(library(tidyverse))

OUT     <- "presentations/smm2026_figs"
CMP_DIR <- Sys.getenv("CMP_DIR", "outputs/joint/compare_s1")
COL_VIS  <- "#D55E00"   # visual (segments), as in make_figures.R
COL_EDNA <- "#0072B2"   # eDNA (stations)

MODES <- c(edna = "eDNA only", visual = "Visual only", joint = "Joint",
           joint_conv = "Joint, c free")
SPECIES <- c(hake = "Hake", humpback = "Humpback", pwsd = "PWSD")

sim <- readRDS(file.path(CMP_DIR, "sim.rds"))
res <- lapply(setNames(names(MODES), names(MODES)), function(m) {
  f <- file.path(CMP_DIR, sprintf("fit_%s.rds", m))
  if (file.exists(f)) readRDS(f)
})
res <- res[!vapply(res, is.null, logical(1))]

src   <- sim$field$locations$source
truth <- sim$field$log_lambda_si
sp    <- sim$field$species
lt    <- sim$meta$lt_species

pts <- map_dfr(names(res), function(m) map_dfr(seq_along(sp), function(s)
  tibble(mode = m, species = sp[s], source = src, truth = truth[, s],
         post = res[[m]]$post_mean[, s]))) |>
  filter(!(mode == "visual" & !species %in% lt)) |>
  mutate(source = factor(source, c("segment", "station"),
                         c("transect segments", "eDNA stations"))) |>
  arrange(source)                                    # stations drawn on top
lab <- map_dfr(res, "field") |>
  filter(!(mode == "visual" & !species %in% lt)) |>
  mutate(label = sprintf("R² %.2f\nCI %.1f", R2, ci_width))

# One row per species, each with square axes on its own range; the visual-only
# panel for a species without visual data is kept empty so columns line up.
modes <- factor(MODES[names(res)], MODES[names(res)])
row_plot <- function(s, first, last) {
  d  <- filter(pts, species == s) |> mutate(mode = factor(MODES[mode], levels(modes)))
  lb <- filter(lab, species == s) |> mutate(mode = factor(MODES[mode], levels(modes)))
  lims <- range(c(d$truth, d$post))
  g <- ggplot(d, aes(truth, post, colour = source)) +
    geom_abline(linetype = 2, colour = "grey45") +
    geom_point(size = 0.6, alpha = 0.45) +
    geom_text(data = lb, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.08, vjust = 1.15, size = 4.6, lineheight = 0.9) +
    facet_wrap(~mode, nrow = 1, drop = FALSE) +
    coord_equal(xlim = lims, ylim = lims) +
    scale_colour_manual(values = c("transect segments" = COL_VIS, "eDNA stations" = COL_EDNA),
                        drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = if (last) "true log density" else NULL, y = SPECIES[[s]], colour = NULL) +
    theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
          axis.title.y = element_text(size = 15, face = "bold"))
  if (!s %in% lt)
    g <- g + geom_text(data = tibble(mode = factor(MODES[["visual"]], levels(modes)),
                                     label = "no visual data"),
                       aes(x = mean(lims), y = mean(lims), label = label),
                       inherit.aes = FALSE, colour = "grey40", size = 5)
  if (!first) g <- g + theme(strip.text = element_blank(), strip.background = element_blank())
  g
}
library(patchwork)
g <- wrap_plots(lapply(seq_along(sp), function(k) row_plot(sp[k], k == 1, k == length(sp))),
                ncol = 1) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT, "results_field.png"), g, width = 12, height = 10, dpi = 150)
if (nzchar(Sys.which("magick")))
  system2("magick", c(shQuote(file.path(OUT, "results_field.png")), "-trim", "+repage",
                      "-bordercolor", "white", "-border", "30",
                      shQuote(file.path(OUT, "results_field.png"))))
cat("Wrote", file.path(OUT, "results_field.png"), "from", CMP_DIR, "\n")
