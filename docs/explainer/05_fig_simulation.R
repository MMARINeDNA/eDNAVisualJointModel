# Figures for the "data simulation" section. Run from the repo root.
source("docs/explainer/00_common.R")
source("R/functions_joint.R")      # also sources functions_visual.R and HSGP4eDNA functions

# ---- one joint simulation: shared field, eDNA stations + transect segments ------------
set.seed(1)
sim <- simulate_joint(seed = 1L, n_stations = 200L, sample_depths = 0)
st  <- sim$edna$design$samples; sg <- sim$visual$design$segments
dom <- LT_DOMAIN

# GP component of the field on a coarse grid (exact GP draw; the model uses an HSGP approximation)
gx <- seq(10, dom$X_km_max - 10, by = 20); gy <- seq(10, dom$Y_km_max - 10, by = 20)
grid <- expand.grid(X = gx, Y = gy); grid$Z_bathy <- seq(100, 3000, length.out = nrow(grid))   # unused: only the GP part is plotted
set.seed(11)
gpf <- simulate_field_bathysp(grid, default_gp_params(), df = 4L)
fld <- bind_cols(grid[c("X", "Y")],
                 as_tibble(gpf$gp_field_si, .name_repair = "minimal") |> setNames(names(default_gp_params()))) |>
  pivot_longer(-c(X, Y), names_to = "species", values_to = "f")
fld$species <- factor(fld$species, names(SP_LAB), SP_LAB)
p1 <- ggplot(fld, aes(X, Y, fill = f)) + geom_raster() + coord_equal() +
  facet_wrap(~species) + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", name = "f(X,Y)") +
  labs(x = "X (km)", y = "Y (km)", title = "Spatial GP component of the latent log-density",
       subtitle = "One anisotropic GP per species: 50 km east-west, 300 km along the coast; marginal SD 2.0 / 1.0 / 1.3") +
  theme_expl()
save_fig(p1, "sim_field_maps", 10, 5)

# ---- bottom-depth effect h_s(Z) and water-column preference ---------------------------
z <- seq(50, 3200, length.out = 200)
bs <- bathy_spline_basis(z, df = 4L, setup = gpf$bathy_spline)$B
h <- map_dfr(names(SP_LAB), function(s) tibble(species = s, Z = z, h = as.numeric(bs %*% BATHY_BETA_TRUE[s, ])))
p2a <- ggplot(h, aes(Z, h, colour = species)) + geom_line(linewidth = 1) +
  scale_colour_manual(values = SP_COL, labels = SP_LAB, name = NULL) +
  labs(x = "Bottom depth Z_bathy (m)", y = "h_s(Z_bathy)  (log density)",
       title = "Bottom-depth effect (natural spline, df 4)") + theme_expl() + theme(legend.position = "bottom")
dp <- c(0, 50, 150, 250, 350, 500)
pref <- tibble(depth = rep(dp, 3), species = rep(names(SP_LAB), each = 6),
               pref = c(0, 1, 3, 2.57, 2.14, 1.5,  0, .67, 2, 1.06, .11, -1.3,  0, -1, -3, -5, -7, -10))
p2b <- ggplot(pref, aes(depth, pref, colour = species)) + geom_line(linewidth = 1) + geom_point() +
  scale_colour_manual(values = SP_COL, labels = SP_LAB, name = NULL) +
  labs(x = "Sampling depth Z_sample (m)", y = "log offset on eDNA copies",
       title = "Water-column preference (known offset)") + theme_expl() + theme(legend.position = "none")
save_fig(p2a + p2b, "sim_depth_effects", 10, 4)

# ---- the joint design: stations and transect segments on the same field ---------------
p3 <- ggplot() +
  geom_segment(data = sg, aes(x = X - seg_l / 2, xend = X + seg_l / 2, y = Y, yend = Y), colour = "grey50", linewidth = 0.4) +
  geom_point(data = st, aes(X, Y), colour = "#b2182b", size = 0.9) +
  coord_equal(xlim = c(0, dom$X_km_max), ylim = c(0, dom$Y_km_max)) +
  labs(x = "X (km)", y = "Y (km)", title = "Simulated joint design",
       subtitle = "200 eDNA stations (red) and 25 east-west\ntransects of 10 km segments (grey), one latent field") +
  theme_expl()
save_fig(p3, "sim_design", 6, 7)

# ---- eDNA observations from the (seed 203-style) depth-preference simulation ----------
es <- readRDS(file.path(EXPL_OUT, "edna_sim.rds"))
es_s <- es$design$samples; ob <- es$observed
qp <- tibble(sample = ob$qpcr_sample_idx, det = ob$qpcr_detect, ct = ob$qpcr_ct) |>
  left_join(tibble(sample = es_s$sample_id, depth = es_s$Z_sample, exp_cp = es$truth$lambda_true_si[, 1] *
                   es$meta$conv_factor[1] * es$truth$zsample_effect[, 1] * es$meta$vol_filtered * es$meta$vol_aliquot / 100), by = "sample")
p4a <- qp |> mutate(bin = cut(log10(exp_cp), seq(-2, 4, 0.5))) |> group_by(bin) |>
  summarise(x = mean(log10(exp_cp)), rate = mean(det), n = n(), .groups = "drop") |>
  ggplot(aes(x, rate)) + geom_point(aes(size = n)) + geom_line() +
  stat_function(fun = function(x) 1 - exp(-es$truth$qpcr_params$kappa * 10^x), colour = "#b2182b") +
  labs(x = "log10 expected copies in the aliquot", y = "qPCR detection rate",
       title = "qPCR hurdle", subtitle = "red: 1 - exp(-kappa x copies), kappa = 0.85") + theme_expl() + theme(legend.position = "none")
p4b <- ggplot(filter(qp, det == 1), aes(log(exp_cp), ct)) + geom_point(alpha = 0.4, size = 0.8) +
  geom_abline(intercept = es$truth$qpcr_params$alpha_ct, slope = -es$truth$qpcr_params$beta_ct, colour = "#b2182b") +
  labs(x = "log expected copies in the reaction", y = "Ct", title = "Ct vs copies",
       subtitle = "red: Ct = alpha - beta ln(copies), beta = 1.44 (about -3.32 per log10)") + theme_expl()
mbp <- as_tibble(ob$mb_reads, .name_repair = "minimal") |> setNames(c(names(SP_LAB), "junk")) |>
  mutate(sample = ob$mb_sample_idx) |> left_join(tibble(sample = es_s$sample_id, depth = es_s$Z_sample), by = "sample") |>
  pivot_longer(names(SP_LAB), names_to = "species", values_to = "reads") |>
  group_by(species, depth) |> summarise(pos = mean(reads > 0), .groups = "drop")
p4c <- ggplot(mbp, aes(depth, pos, colour = species)) + geom_line(linewidth = 1) + geom_point() +
  scale_colour_manual(values = SP_COL, labels = c(hake = "hake", humpback = "humpback", pwsd = "PWSD"), name = NULL) +
  labs(x = "Sampling depth (m)", y = "share of MB replicates with reads",
       title = "Metabarcoding positives by depth", subtitle = "the depth preferences are visible in the raw data") +
  theme_expl() + theme(legend.position = "bottom")
save_fig(p4a + p4b + p4c, "sim_edna_obs", 14, 4.4)

# ---- visual observations: detection function, distances, group sizes -------------------
sp <- lt_species_params()
dd <- bind_rows(lapply(c("humpback", "pwsd"), function(s) transform(sim$visual$observed[[s]]$obs, species = s)))
xg <- seq(0, 5, length.out = 100)
mix_g <- function(s) {   # density of detected distances: mean over group sizes of g(x | s)
  p <- sp[[s]]
  g <- if (p$use_size_covar == 1L)
    rowMeans(sapply(sample(p$group_size_pool, 400, TRUE), function(k)
      exp(-xg^2 / (2 * (p$sigma_det * exp(p$beta_size_truth * (k - p$s_centre)))^2)))) else
    exp(-xg^2 / (2 * p$sigma_det^2))
  tibble(species = s, x = xg, dens = g / (sum(g) * (xg[2] - xg[1])))
}
det_curve <- bind_rows(mix_g("humpback"), mix_g("pwsd"))
p5a <- ggplot(dd, aes(distance, after_stat(density))) + geom_histogram(binwidth = 0.5, boundary = 0, closed = "left", fill = "grey80", colour = "white") +
  geom_line(data = det_curve, aes(x, dens), colour = "#b2182b", linewidth = 1) +
  facet_wrap(~species, scales = "free_y") + labs(x = "perpendicular distance (km)", y = "density",
  title = "Detected distances against the true detection function",
  subtitle = "red: truth, averaged over group sizes (half-normal; PWSD scale grows with group size)") + theme_expl()
pools <- bind_rows(lapply(c("humpback", "pwsd"), function(s) tibble(species = s, size = sp[[s]]$group_size_pool)))
p5b <- ggplot(pools, aes(size)) + geom_histogram(bins = 30, fill = "grey60", colour = "white") +
  scale_x_log10() + facet_wrap(~species, scales = "free") +
  labs(x = "group size (log scale)", y = "sightings in the pool", title = "Empirical group-size pools (data/grpsz)") + theme_expl()
save_fig(p5a / p5b, "sim_visual_obs", 9, 7)
cat("done\n")
