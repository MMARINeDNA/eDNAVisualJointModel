# Diagram: how sea state (Beaufort) is simulated and how it changes detectability.
source("docs/explainer/00_common.R")
source("R/functions_visual.R")      # LT_DOMAIN, BFT_CENTRE and draw_sea_state() as the simulator uses them
suppressMessages(library(patchwork))
DOM <- list(X = LT_DOMAIN$X_km_max, Y = LT_DOMAIN$Y_km_max)
beta <- c(humpback = -0.08, pwsd = -0.15); sig <- c(humpback = 2.5, pwsd = 1.5)

# A. flow of the simulation
boxes <- tibble(x = c(1, 3.3, 5.6, 7.9), label = c(
  "Segment sea state\nBeaufort b\n3 + north-south trend\n+ east-west trend + noise\n(clipped to 0 to 6)",
  "Detection scale\nsigma = sigma_det\n x exp(beta_size (s - s_bar))\n x exp(beta_bft (b - 3))",
  "Detection probability\ng(x) = exp(-x^2 / 2 sigma^2)\nat perpendicular\ndistance x",
  "Group is kept\nwith probability g(x)\n(density itself does\nnot depend on b)"))
pA <- ggplot(boxes) +
  geom_rect(aes(xmin = x - 1, xmax = x + 1, ymin = 0, ymax = 1), fill = "grey95", colour = "grey40") +
  geom_text(aes(x, 0.5, label = label), size = 3.2, lineheight = 1.05) +
  annotate("segment", x = c(2.05, 4.35, 6.65), xend = c(2.25, 4.55, 6.85), y = 0.5, yend = 0.5,
           arrow = arrow(length = unit(0.12, "inches")), linewidth = 0.8) +
  coord_cartesian(xlim = c(-0.1, 9.1), ylim = c(-0.05, 1.05)) + theme_void() +
  labs(title = "Sea state and detectability in the simulator")

# B. the Beaufort field (expected value, before the N(0, 0.5) noise)
g <- expand.grid(X = seq(0, DOM$X, length.out = 60), Y = seq(0, DOM$Y, length.out = 120))
# expected value of draw_sea_state(): BFT_CENTRE + 1.0 (2Y/Ymax - 1) + 0.75 (1 - 2X/Xmax), clipped to [0, 6]
g$b <- pmin(pmax(BFT_CENTRE + 1.0 * (2 * g$Y / DOM$Y - 1) + 0.75 * (1 - 2 * g$X / DOM$X), 0), 6)
pB <- ggplot(g, aes(X, Y, fill = b)) + geom_raster() + coord_equal() +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Beaufort") +
  labs(x = "X (km)", y = "Y (km)", title = "Sea state across the survey",
       subtitle = "expected value; segments add N(0, 0.5) noise") + theme_expl()

# C. detection functions at Beaufort 1, 3, 5
x <- seq(0, 5, length.out = 200)
cur <- expand_grid(species = names(beta), bft = c(1, 3, 5), x = x) |>
  mutate(sigma = sig[species] * exp(beta[species] * (bft - BFT_CENTRE)), g = exp(-x^2 / (2 * sigma^2)),
         species = factor(species, names(beta), SP_LAB[names(beta)]))
lab <- cur |> distinct(species, bft, sigma)
pC <- ggplot(cur, aes(x, g, colour = factor(bft))) + geom_line(linewidth = 1) + facet_wrap(~species) +
  scale_colour_viridis_d(option = "magma", end = 0.8, direction = -1, name = "Beaufort") +
  geom_text(data = lab, aes(x = 4.9, y = 0.97 - 0.1 * match(bft, c(1, 3, 5)),
                            label = sprintf("Beaufort %d: sigma = %.2f km", bft, sigma)), hjust = 1, size = 3, colour = "black", show.legend = FALSE) +
  labs(x = "perpendicular distance (km)", y = "detection probability g(x)",
       title = "Detection function by sea state",
       subtitle = "PWSD at the reference group size\n(larger groups widen sigma further)") + theme_expl()

# D. average detection probability in a 5 km strip against Beaufort
bb <- seq(0, 6, 0.25)
pd <- expand_grid(species = names(beta), b = bb) |>
  mutate(sigma = sig[species] * exp(beta[species] * (b - BFT_CENTRE)),
         p = sigma * sqrt(pi / 2) * (2 * pnorm(5 / sigma) - 1) / 5,
         species = factor(species, names(beta), SP_LAB[names(beta)]))
pD <- ggplot(pd, aes(b, p, colour = species)) + geom_line(linewidth = 1) +
  scale_colour_manual(values = c(setNames(SP_COL[c("humpback", "pwsd")], SP_LAB[c("humpback", "pwsd")])), name = NULL) +
  labs(x = "Beaufort", y = "mean detection probability in a 5 km strip",
       title = "Average detection probability",
       subtitle = "beta = -0.08 (humpback), -0.15 (PWSD)\nper Beaufort unit") + theme_expl() +
  theme(legend.position = "bottom")

p <- pA / (pB + pC + pD + plot_layout(widths = c(0.8, 1.7, 1))) + plot_layout(heights = c(0.55, 1.5))
save_fig(p, "sim_sea_state", 15, 8)
print(pd |> filter(b %in% c(1, 3, 5)) |> mutate(p = round(p, 3)) |> select(species, b, sigma, p))
