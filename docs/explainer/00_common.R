# Shared setup for the explainer scripts. Run from the repo root.
suppressMessages({ library(cmdstanr); library(posterior); library(tidyverse); library(patchwork) })
HSGP4EDNA_DIR <- "external/HSGP4eDNA"
EXPL_OUT <- "outputs/explainer"           # fits and intermediates (git-ignored)
FIG_DIR  <- "docs/figures"
dir.create(EXPL_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

SP_COL <- c(hake = "#1b7837", humpback = "#2166ac", pwsd = "#b2182b")
SP_LAB <- c(hake = "Pacific hake", humpback = "Humpback whale", pwsd = "Pacific white-sided dolphin")
theme_expl <- function() theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
save_fig <- function(p, name, w = 9, h = 5) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 130)
  invisible(p)
}
