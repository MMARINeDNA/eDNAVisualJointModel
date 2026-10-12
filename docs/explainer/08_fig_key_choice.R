# Real-data distances with half-normal and hazard-rate fits (Distance package, no covariates).
source("docs/explainer/00_common.R")
# Maximum-likelihood fits of the truncated detection density f(x) = g(x) / ESW, no covariates,
# written out directly (the same density as stan/hsgp_visual.stan) rather than via a package.
vis <- readRDS("data/processed/muri_visual.rds")
g_hn <- function(x, par) exp(-x^2 / (2 * exp(par[1])^2))
g_hr <- function(x, par) -expm1(-(exp(par[1]) / pmax(x, 1e-6))^exp(par[2]))
nll <- function(par, x, w, gfun) {
  xx <- (seq_len(2000) - 0.5) / 2000 * w
  -sum(log(gfun(x, par))) + length(x) * log(mean(gfun(xx, par)) * w)
}
rows <- list(); curves <- list(); aic <- list()
for (sp in c("humpback", "pwsd")) {
  w <- vis$meta$w[[sp]]
  d <- vis$observed[[sp]]$obs; d <- d[d$distance <= w, ]
  rows[[sp]] <- tibble(species = sp, distance = d$distance)
  xg <- seq(0, w, length.out = 100)
  for (key in c("hn", "hr")) {
    g <- if (key == "hn") g_hn else g_hr
    start <- if (key == "hn") log(2) else c(log(2), log(2))
    o <- if (key == "hn") {
      r <- optimize(function(p) nll(p, d$distance, w, g), c(-3, 4)); list(par = r$minimum, value = r$objective)
    } else optim(start, nll, x = d$distance, w = w, gfun = g, method = "Nelder-Mead", control = list(maxit = 2000))
    gg <- g(xg, o$par)
    avgp <- mean(g((seq_len(4000) - 0.5) / 4000 * w, o$par))
    curves[[paste(sp, key)]] <- tibble(species = sp, key = c(hn = "half-normal", hr = "hazard-rate")[key], x = xg,
                                       dens = gg / (sum(gg) * (xg[2] - xg[1])))
    aic[[paste(sp, key)]] <- tibble(species = sp, key = key, n = nrow(d), AIC = 2 * o$value + 2 * length(o$par),
                                    avg_p = avgp, scale_km = exp(o$par[1]), shape = if (key == "hr") exp(o$par[2]) else NA)
  }
}
dd <- bind_rows(rows); cv <- bind_rows(curves); ac <- bind_rows(aic)
ac$delta <- ac$AIC - ave(ac$AIC, ac$species, FUN = min)
saveRDS(ac, file.path(EXPL_OUT, "visual_key_choice.rds"))
p <- ggplot(dd, aes(distance, after_stat(density))) + geom_histogram(binwidth = 0.5, fill = "grey80", colour = "white") +
  geom_line(data = cv, aes(x, dens, colour = key), linewidth = 1) +
  scale_colour_manual(values = c("half-normal" = "#2166ac", "hazard-rate" = "#b2182b"), name = NULL) +
  facet_wrap(~species, scales = "free") + labs(x = "perpendicular distance (km)", y = "density",
  title = "Real data: detected distances with half-normal and hazard-rate fits",
  subtitle = "no covariates; truncated at the 95th percentile per species") + theme_expl() + theme(legend.position = "bottom")
save_fig(p, "visual_key_choice", 10, 4.2)
print(ac)
