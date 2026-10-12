# Figures for the visual model: simulated survey (with truth) and real MURI data.
# Needs outputs/explainer/visual_{sim,real}_<species>_{standata,fit}.rds.
source("docs/explainer/00_common.R")
source("R/functions_visual.R")

hn_g <- function(x, sig) exp(-x^2 / (2 * sig^2))
hr_g <- function(x, sig, b) -expm1(-(sig / pmax(x, 1e-6))^b)
qig <- function(p, shape, scale) 1 / qgamma(1 - p, shape, rate = scale)

visual_figs <- function(tag, sp, sd, fit, seg, obs_all = NULL, truth = NULL, w) {
  draws <- as_draws_matrix(fit$draws())
  v <- function(name) as.numeric(draws[, name])
  cols <- function(prefix) which(startsWith(colnames(draws), prefix))
  nseg <- sd$n_seg
  lab <- SP_LAB[[sp]]; out <- list()
  # ---- diagnostics
  dg <- fit$diagnostic_summary(quiet = TRUE)
  keyv <- c("mu_sp", "log_sigma", "beta_size", "gp_sigma[1]", "gp_l[1]", "gp_l[2]",
            if (sd$det_key == 1L) "b_shape", if (sd$use_bft_covar == 1L) "beta_bft[1]", "mean_group_size", "D_mean")
  keyv <- intersect(keyv, colnames(draws))
  sm <- summarise_draws(draws[, keyv], mean, ~quantile(.x, c(.05, .5, .95)), rhat, ess_bulk)
  out$diag <- list(div = sum(dg$num_divergent), td = sum(dg$num_max_treedepth), ebfmi = dg$ebfmi,
                   max_rhat = max(sm$rhat, na.rm = TRUE), min_ess = min(sm$ess_bulk, na.rm = TRUE),
                   runtime_min = max(fit$time()$chains$total) / 60, M = sd$M, n = sd$n, n_seg = nseg, summary = sm)

  # ---- field at segments
  ll <- draws[, cols("log_lambda["), drop = FALSE]
  ll <- ll[, !startsWith(colnames(ll), "log_lambda_pred"), drop = FALSE]
  pm <- colMeans(ll)
  seg$post_log_lambda <- pm
  if (!is.null(truth)) {
    seg$true_log_lambda <- truth$log_lambda
    out$r2 <- cor(pm, truth$log_lambda)^2
    p <- ggplot(seg, aes(true_log_lambda, post_log_lambda)) + geom_point(alpha = 0.35, size = 0.8) +
      geom_abline(colour = "#b2182b") + labs(x = "true log density (animals/km²)",
        y = "posterior mean log density", title = sprintf("%s: latent field at segments", lab),
        subtitle = sprintf("R² = %.2f", out$r2)) + theme_expl()
    save_fig(p, sprintf("visual_%s_%s_field", tag, sp), 5, 4.5)
  }
  # map of posterior mean density with sightings
  o <- sd$obs_df
  pmap <- ggplot(seg, aes(X, Y, colour = exp(post_log_lambda))) + geom_point(size = 0.6) +
    scale_colour_viridis_c(trans = "log10", name = "animals / km²") +
    { if (!is.null(o)) geom_point(data = o, aes(X, Y), inherit.aes = FALSE, shape = 1, colour = "#b2182b", size = 1.6) } +
    labs(x = "X (km)", y = "Y (km)", title = sprintf("%s: posterior mean density along the effort", lab),
         subtitle = "red circles: detections within the truncation distance") + coord_equal() + theme_expl()
  out$map <- pmap

  # ---- posterior predictive: counts per transect
  pp <- draws[, cols("pp_seg_count["), drop = FALSE]
  tr <- if (n_distinct(seg$transect_id) > 100) as.integer(cut(seg$Y, seq(0, max(seg$Y) + 100, by = 100))) else seg$transect_id
  agg <- map_dfr(sort(unique(tr)), function(t) {
    j <- which(tr == t); s <- rowSums(pp[, j, drop = FALSE])
    tibble(transect = t, obs = sum(sd$seg_count[j]), lo = quantile(s, .05), med = median(s), hi = quantile(s, .95),
           Y = mean(seg$Y[j]))
  })
  out$count_unit <- if (n_distinct(seg$transect_id) > 100) "100 km latitude band" else "transect line"
  out$ppc_counts <- agg
  # ---- posterior predictive: distances (density of detected distances)
  xg <- seq(0, w, length.out = 120)
  nd <- min(150, nrow(draws)); di <- round(seq(1, nrow(draws), length.out = nd))
  pkc <- cols("pk[")
  curves <- sapply(di, function(i) {
    pk <- as.numeric(draws[i, pkc]); K <- length(pk)
    ls <- v("log_sigma")[i] + (if (sd$use_size_covar == 1L) v("beta_size")[i] * (seq_len(K) - sd$s_centre) else 0)
    if (sd$use_bft_covar == 1L) ls <- ls + v("beta_bft[1]")[i] * mean(sd$bft_level[sd$seg_bft])
    ks <- if (sd$use_size_covar == 1L) which(pk > 1e-4) else 1
    dens <- rowSums(sapply(ks, function(k) (if (sd$det_key == 1L) hr_g(xg, exp(ls[k]), exp(v("log_b[1]")[i])) else
      hn_g(xg, exp(ls[k]))) * (if (sd$use_size_covar == 1L) pk[k] else 1)))
    dens / (sum(dens) * (xg[2] - xg[1]))
  })
  out$ppc_dist <- list(x = xg, lo = apply(curves, 1, quantile, .05), med = apply(curves, 1, median),
                       hi = apply(curves, 1, quantile, .95), obs = sd$x)
  # ---- hyperparameters vs prior
  cs <- sd$coord_scale
  hp <- tibble(param = c("lx (km)", "ly (km)", "gp_sigma", "mu_sp"),
    var = c("gp_l[1]", "gp_l[2]", "gp_sigma[1]", "mu_sp"),
    plo = c(qig(.025, sd$prior_gp_l_shape[1], sd$prior_gp_l_scale[1]) * cs[1], qig(.025, sd$prior_gp_l_shape[2], sd$prior_gp_l_scale[2]) * cs[2],
            qgamma(.025, sd$prior_gp_sigma_shape, sd$prior_gp_sigma_rate), qnorm(.025, sd$prior_mu_sp_mu, sd$prior_mu_sp_sig)),
    phi = c(qig(.975, sd$prior_gp_l_shape[1], sd$prior_gp_l_scale[1]) * cs[1], qig(.975, sd$prior_gp_l_shape[2], sd$prior_gp_l_scale[2]) * cs[2],
            qgamma(.975, sd$prior_gp_sigma_shape, sd$prior_gp_sigma_rate), qnorm(.975, sd$prior_mu_sp_mu, sd$prior_mu_sp_sig))) |>
    mutate(med = map_dbl(var, ~ median(v(.x))), lo = map_dbl(var, ~ quantile(v(.x), .025)), hi = map_dbl(var, ~ quantile(v(.x), .975)))
  if (!is.null(truth$hyper)) hp$truth <- truth$hyper
  out$hyper <- hp
  gpp <- attr(sd, "gp_prior")
  if (!is.null(gpp)) out$floor <- floor_check(cbind(v("gp_l[1]"), v("gp_l[2]")), gpp)
  # ---- prior predictive: expected number of detected groups
  set.seed(5)
  coords_km <- as.matrix(seg[, c("X", "Y")]); B <- sd$B_bathy
  idx <- if (nseg > 600) sort(sample(nseg, 600)) else seq_len(nseg)   # thin for the exact GP
  sc <- nseg / length(idx)
  pg <- replicate(150, {
    gsig <- rgamma(1, sd$prior_gp_sigma_shape, sd$prior_gp_sigma_rate)
    ell <- c(1 / rgamma(1, sd$prior_gp_l_shape[1], rate = sd$prior_gp_l_scale[1]),
             1 / rgamma(1, sd$prior_gp_l_shape[2], rate = sd$prior_gp_l_scale[2])) * cs
    f <- as.numeric(MASS::mvrnorm(1, rep(0, length(idx)), aniso_cov_2d(coords_km[idx, ], gsig, ell[1], ell[2])))
    mu <- rnorm(1, sd$prior_mu_sp_mu, sd$prior_mu_sp_sig)
    h <- as.numeric(B[idx, , drop = FALSE] %*% rnorm(sd$K_bathy, 0, sd$prior_beta_bathy_sig))
    ls <- rnorm(1, sd$log_sigma_prior_mean, sd$log_sigma_prior_sd)
    esw <- if (sd$det_key == 1L) mean(hr_g(seq(0.01, w, length.out = 300), exp(min(ls, 8)), exp(rnorm(1, sd$log_b_prior_mean, sd$log_b_prior_sd)))) * w else
      exp(ls) * sqrt(pi / 2) * erf_(w / (sqrt(2) * exp(ls)))
    ms <- sd$mu_log_prior_mean   # mean group size ~ exp(mu_log) for the lognormal species; NB mean otherwise
    msz <- if (sd$group_size_dist == 1L) exp(rnorm(1, sd$mu_log_prior_mean, sd$mu_log_prior_sd)) else
      rgamma(1, sd$mu_s_prior_shape, sd$mu_s_prior_rate)
    sum(exp(mu + f + h)) * sc / msz * 2 * mean(sd$seg_l) * esw
  })
  out$prior_total <- list(draws = pg, obs = sum(sd$seg_count))
  out$post_total <- rowSums(pp)
  out
}
erf_ <- function(x) 2 * pnorm(x * sqrt(2)) - 1
aniso_cov_2d <- function(coords, sigma, lx, ly, jitter = 1e-6) {
  dx <- outer(coords[, 1], coords[, 1], "-") / lx; dy <- outer(coords[, 2], coords[, 2], "-") / ly
  sigma^2 * exp(-0.5 * (dx^2 + dy^2)) + diag(jitter, nrow(coords))
}

plot_visual_set <- function(res, tag, sp) {
  lab <- SP_LAB[[sp]]
  p1 <- ggplot(res$ppc_counts, aes(Y)) + geom_pointrange(aes(y = med, ymin = lo, ymax = hi), colour = "#2166ac", size = 0.3) +
    geom_point(aes(y = obs), colour = "#b2182b", shape = 4, size = 2.2, stroke = 1.1) +
    labs(x = "position Y (km)", y = paste("detected groups per", res$count_unit),
         title = sprintf("%s: posterior predictive counts", lab), subtitle = "blue: replicated (90%); red x: observed") + theme_expl()
  d <- res$ppc_dist
  p2 <- ggplot() + geom_histogram(data = tibble(x = d$obs), aes(x, after_stat(density)), binwidth = max(d$x) / 10,
          boundary = 0, closed = "left", fill = "grey80", colour = "white") +
    geom_ribbon(data = tibble(x = d$x, lo = d$lo, hi = d$hi), aes(x, ymin = lo, ymax = hi), fill = "#2166ac", alpha = 0.3) +
    geom_line(data = tibble(x = d$x, m = d$med), aes(x, m), colour = "#2166ac", linewidth = 1) +
    labs(x = "perpendicular distance (km)", y = "density", title = sprintf("%s: detection function", lab),
         subtitle = "blue: fitted density of detected distances (90%)") + theme_expl()
  pt <- tibble(total = c(res$prior_total$draws, res$post_total), src = rep(c("prior", "posterior"), c(length(res$prior_total$draws), length(res$post_total))))
  p3 <- ggplot(pt, aes(total, fill = src)) + geom_density(aes(y = after_stat(scaled)), alpha = 0.6, adjust = 1.2) +
    geom_vline(xintercept = res$prior_total$obs, colour = "#b2182b", linewidth = 1) + scale_x_log10() +
    scale_fill_manual(values = c(prior = "grey55", posterior = "#2166ac"), name = NULL) +
    labs(x = "total detected groups (log scale)", y = "density (scaled)", title = "Prior vs posterior predictive total",
         subtitle = "red: observed") + theme_expl() + theme(legend.position = "bottom")
  save_fig(p1 + p2 + p3, sprintf("visual_%s_%s_ppc", tag, sp), 14, 4.4)
  h <- res$hyper |> mutate(param = factor(param, unique(param)))
  p4 <- ggplot(h, aes(x = 1)) + geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = plo, ymax = phi), fill = "grey88", inherit.aes = FALSE) +
    geom_pointrange(aes(y = med, ymin = lo, ymax = hi), colour = "#2166ac") + facet_wrap(~param, scales = "free_y", nrow = 1) +
    { if ("truth" %in% names(h)) geom_point(aes(y = truth), colour = "#b2182b", shape = 4, size = 3, stroke = 1.3) } +
    labs(x = NULL, y = NULL, title = sprintf("%s: GP hyperparameters", lab),
         subtitle = paste("blue: posterior 95%; grey: central 95% of the prior", if ("truth" %in% names(h)) "; red x: truth" else "")) +
    theme_expl() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  save_fig(p4, sprintf("visual_%s_%s_hyper", tag, sp), 10, 3.2)
  save_fig(res$map, sprintf("visual_%s_%s_map", tag, sp), 6, 7.5)
}

# ---- simulated survey ------------------------------------------------------------------------
sim_f <- file.path(EXPL_OUT, "visual_sim.rds")
if (file.exists(sim_f)) {
  sim <- readRDS(sim_f); res_sim <- list()
  for (sp in c("humpback", "pwsd")) {
    ff <- file.path(EXPL_OUT, sprintf("visual_sim_%s_fit.rds", sp))
    if (!file.exists(ff)) next
    sd <- readRDS(file.path(EXPL_OUT, sprintf("visual_sim_%s_standata.rds", sp))); fit <- readRDS(ff)
    seg <- sim$design$segments
    o <- sim$observed[[sp]]$obs; o <- o[o$distance <= sim$meta$w, ]
    sd$obs_df <- tibble(X = seg$X[o$seg_id], Y = seg$Y[o$seg_id])
    r <- visual_figs("sim", sp, sd, fit, seg, truth = list(log_lambda = sim$truth$log_lambda[, sp],
                   hyper = c(sim$truth$gp_params[[sp]]$lx, sim$truth$gp_params[[sp]]$ly, sim$truth$gp_params[[sp]]$sigma,
                             sim$truth$gp_params[[sp]]$mu)), w = sim$meta$w)
    plot_visual_set(r, "sim", sp)
    # parameter recovery table
    p <- sim$truth$sp_params[[sp]]; g <- sim$truth$gp_params[[sp]]
    draws <- as_draws_matrix(fit$draws())
    tr <- tibble(param = c("sigma_c (km)", "D mean (animals/km2)", "lx (km)", "ly (km)", "gp_sigma", "mean group size"),
                 var = c("sigma_c", "D_mean", "gp_l[1]", "gp_l[2]", "gp_sigma[1]", "mean_group_size"),
                 truth = c(p$sigma_det, mean(exp(sim$truth$log_lambda[, sp])), g$lx, g$ly, g$sigma, p$mean_group_size))
    tr$med <- map_dbl(tr$var, ~ median(draws[, .x])); tr$lo <- map_dbl(tr$var, ~ quantile(draws[, .x], .025))
    tr$hi <- map_dbl(tr$var, ~ quantile(draws[, .x], .975)); tr$covered <- tr$truth >= tr$lo & tr$truth <= tr$hi
    r$recovery <- tr; res_sim[[sp]] <- r[c("diag", "r2", "recovery", "floor")]
  }
  saveRDS(res_sim, file.path(EXPL_OUT, "visual_sim_results.rds"))
}

# ---- real data (default budget, and the refit at twice the budget: real800 before v1.4, real1600 after) ---------------------------------
vis_f <- "data/processed/muri_visual.rds"
res_real <- list()
if (file.exists(vis_f)) {
  vis <- readRDS(vis_f)
  for (tag in c("real", "real800", "real1600")) for (sp in c("humpback", "pwsd")) {
    ff <- file.path(EXPL_OUT, sprintf("visual_%s_%s_fit.rds", tag, sp))
    if (!file.exists(ff)) next
    sd <- readRDS(file.path(EXPL_OUT, sprintf("visual_%s_%s_standata.rds", tag, sp))); fit <- readRDS(ff)
    seg <- vis$design$segments
    o <- vis$observed[[sp]]$obs; o <- o[o$distance <= vis$meta$w[[sp]], ]
    sd$obs_df <- tibble(X = seg$X[o$seg_id], Y = seg$Y[o$seg_id])
    r <- visual_figs(tag, sp, sd, fit, seg, w = vis$meta$w[[sp]])
    plot_visual_set(r, tag, sp)
    res_real[[paste(tag, sp, sep = "_")]] <- r[c("diag", "floor", "hyper")]
  }
  saveRDS(res_real, file.path(EXPL_OUT, "visual_real_results.rds"))
}
cat("done\n")
